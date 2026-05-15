import os
import argparse
import pandas as pd
import gzip  
from collections import defaultdict


# ==========================================
# 1. PARSER MODULE (VCF Extraction)
# ==========================================
class VCFExtractor:
    def __init__(self, pgx_variants_path):
        self.variants_df = pd.read_csv(pgx_variants_path, sep='\t')
        # Создаем словарь для быстрого поиска: (chrom, pos) -> dict
        self.target_coords = {}
        for _, row in self.variants_df.iterrows():
            self.target_coords[(str(row['chrom']).replace('chr', ''), str(row['pos']))] = row.to_dict()
            self.target_coords[(str(row['chrom']), str(row['pos']))] = row.to_dict()

    def parse_vcf(self, vcf_path):
        extracted = []
        samples = []
        
        # Функция для умного открытия файлов (сжатых и несжатых)
        def smart_open(filepath):
            with open(filepath, 'rb') as test_f:
                magic_number = test_f.read(2)
            if magic_number == b'\x1f\x8b':
                return gzip.open(filepath, 'rt', encoding='utf-8')
            else:
                return open(filepath, 'r', encoding='utf-8')

        # Fallback parser (pure python)
        with smart_open(vcf_path) as f:
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 5: 
                    continue
                    
                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
                
                if (chrom, pos) in self.target_coords:
                    var_info = self.target_coords[(chrom, pos)]
                    
                    for i, sample in enumerate(samples):
                        if 9 + i < len(parts):
                            gt_field = parts[9 + i].split(':')[0]
                        else:
                            gt_field = '.'
                            
                        gt_status = self._map_gt(gt_field)
                        
                        extracted.append({
                            'sample': sample,
                            'gene': var_info['gene'],
                            'marker_id': var_info['marker_id'],
                            'marker_type': var_info['marker_type'],
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'raw_gt': gt_field,
                            'genotype': gt_status
                        })
        return pd.DataFrame(extracted), samples

    def _map_gt(self, gt):
        if gt in ['0/0', '0|0']: return 'HOM_REF'
        elif gt in ['0/1', '1/0', '0|1', '1|0']: return 'HET'
        elif gt in ['1/1', '1|1']: return 'HOM_ALT'
        elif '.' in gt: return 'NO_CALL'
        else: return 'AMBIGUOUS'

# 2. CALLER MODULE (Diplotype & Genotype)
# ==========================================
class PGxCaller:
    def __init__(self, star_alleles_path):
        self.star_alleles = pd.read_csv(star_alleles_path, sep='\t')
        
    def call_genes(self, extracted_df):
        calls = []
        coverage_stats = []
        
        for sample, sample_df in extracted_df.groupby('sample'):
            # Считаем кавередж
            total_markers = len(sample_df)
            called_markers = len(sample_df[sample_df['genotype'] != 'NO_CALL'])
            coverage_stats.append({
                'sample': sample,
                'total_targeted_markers': total_markers,
                'successfully_called': called_markers,
                'call_rate': called_markers / total_markers if total_markers > 0 else 0
            })
            
            for gene, gene_df in sample_df.groupby('gene'):
                marker_type = gene_df['marker_type'].iloc[0]
                
                if marker_type.lower() == 'snp':
                    # Для SNP просто передаем генотип
                    for _, row in gene_df.iterrows():
                        calls.append({
                            'sample': sample,
                            'gene': gene,
                            'marker_id': row['marker_id'],
                            'call_type': 'SNP',
                            'result': row['genotype'], # HOM_REF, HET, HOM_ALT
                            'confidence': 'HIGH' if row['genotype'] != 'NO_CALL' else 'LOW',
                            'warnings': '' if row['genotype'] != 'NO_CALL' else 'Missing data'
                        })
                else:
                    # Логика Star Alleles
                    diplotype, confidence, warning = self._call_star_allele(gene, gene_df)
                    calls.append({
                        'sample': sample,
                        'gene': gene,
                        'marker_id': 'ALL',
                        'call_type': 'STAR',
                        'result': diplotype,
                        'confidence': confidence,
                        'warnings': warning
                    })
                    
        return pd.DataFrame(calls), pd.DataFrame(coverage_stats)

    def _call_star_allele(self, gene, gene_df):
        # Базовая эвристика вызова star-аллелей (без фазирования)
        # Получаем мутации пациента (HET или HOM_ALT)
        mutated_markers = gene_df[gene_df['genotype'].isin(['HET', 'HOM_ALT'])]['marker_id'].tolist()
        hom_alt_markers = gene_df[gene_df['genotype'] == 'HOM_ALT']['marker_id'].tolist()
        
        gene_stars = self.star_alleles[self.star_alleles['gene'] == gene]
        if gene_stars.empty:
            return '*1/*1', 'LOW', 'No reference star alleles found'

        detected_stars = set()
        
        for star, star_df in gene_stars.groupby('star'):
            required = set(star_df['marker_id'].tolist())
            if required and required.issubset(set(mutated_markers)):
                detected_stars.add(star)
                
        # Формируем дилотипы (очень упрощенно)
        detected_stars = list(detected_stars)
        if len(detected_stars) == 0:
            return '*1/*1', 'HIGH', ''
        elif len(detected_stars) == 1:
            # Проверяем, это гетерозигота или гомозигота
            star = detected_stars[0]
            req_markers = set(gene_stars[gene_stars['star'] == star]['marker_id'])
            if req_markers.issubset(set(hom_alt_markers)):
                return f"{star}/{star}", 'MODERATE', 'Inferred homozygous from unphased data'
            else:
                return f"*1/{star}", 'MODERATE', 'Inferred heterozygous'
        else:
            # Два или более. Берем первые два (по-хорошему нужен phasing)
            return f"{detected_stars[0]}/{detected_stars[1]}", 'LOW', 'Unphased multiple star alleles detected'

# ==========================================
# 3. RULES ENGINE MODULE
# ==========================================
class RulesEngine:
    def __init__(self, rules_dir):
        self.rules_snp = pd.read_csv(os.path.join(rules_dir, 'rulessnp_rules.tsv'), sep='\t')
        self.rules_star = pd.read_csv(os.path.join(rules_dir, 'rulesstar_alleles.tsv'), sep='\t')
        self.rules_comb_ssnp = pd.read_csv(os.path.join(rules_dir, 'rules_comb_star_alleles_snp.tsv'), sep='\t')
        self.rules_comb_ss = pd.read_csv(os.path.join(rules_dir, 'rules_comb_star_alleles.tsv'), sep='\t')

        # Создаем словарь для автоматического перевода Диплотип -> Фенотип
        # Ключ: (gene, diplotype), Значение: phenotype
        self.phenotype_map = {}
        for _, row in self.rules_star.iterrows():
            if pd.notna(row.get('phenotype')):
                self.phenotype_map[(row['gene'], row['diplotype'])] = row['phenotype']

    def apply_rules(self, calls_df):
        reports = []
        
        for sample, sample_calls in calls_df.groupby('sample'):
            # Кэшируем вызовы для быстрого доступа
            snp_calls = sample_calls[sample_calls['call_type'] == 'SNP'].set_index('marker_id')['result'].to_dict()
            star_calls = sample_calls[sample_calls['call_type'] == 'STAR'].set_index('gene')['result'].to_dict()
            
            # ВЫЧИСЛЯЕМ ФЕНОТИПЫ пациента на основе его диплотипов
            star_phenotypes = {}
            for gene, diplotype in star_calls.items():
                star_phenotypes[gene] = self.phenotype_map.get((gene, diplotype), 'Unknown')
                
            processed_drugs = set()

            # --- ЭТАП 1: Комбинированные правила (Самый высокий приоритет) ---
            
            # 1.1 Star + Star (Здесь сопоставляем ФЕНОТИПЫ, т.к. в файле phenotype1 и phenotype2)
            for _, rule in self.rules_comb_ss.iterrows():
                drug = rule['drug']
                g1, g2 = rule['gene1'], rule['gene2']
                
                if drug not in processed_drugs:
                    if (g1 in star_phenotypes and star_phenotypes[g1] == rule['phenotype1'] and 
                        g2 in star_phenotypes and star_phenotypes[g2] == rule['phenotype2']):
                        
                        result_str = f"{star_calls[g1]} ({star_phenotypes[g1]}) + {star_calls[g2]} ({star_phenotypes[g2]})"
                        phenotype_str = f"{star_phenotypes[g1]}+{star_phenotypes[g2]}"
                        reports.append(self._build_report_row(sample, rule, drug, 'COMBINED_STAR_STAR', f"{g1}+{g2}", result_str, phenotype_str))
                        processed_drugs.add(drug)

	# 1.2 Star + SNP (Здесь сопоставляем ДИПЛОТИПЫ, т.к. в файле diplotype)
            for _, rule in self.rules_comb_ssnp.iterrows():
                # .strip() убивает все случайные невидимые пробелы из TSV-таблиц
                drug = str(rule['drug']).strip()
                snp_marker = str(rule.get('marker_id', '')).strip()
                star_gene = str(rule.get('star_gene', '')).strip()
                req_gt = str(rule.get('genotype', '')).strip()
                req_dip = str(rule.get('diplotype', '')).strip()
                
                if drug not in processed_drugs:
                    # Безопасно достаем статусы пациента (если маркера нет, будет пустая строка '')
                    pat_gt = str(snp_calls.get(snp_marker, '')).strip()
                    pat_dip = str(star_calls.get(star_gene, '')).strip()

                    # СПЕЦИАЛЬНЫЙ ДЕБАГ: Выводим в консоль то, что видит пайплайн
                    if drug.lower() == 'warfarin' and sample == 'TEST_PATIENT_COMB':
                        print(f"\n--- DEBUG {sample} ---")
                        print(f"Правило требует: {snp_marker} = '{req_gt}' И {star_gene} = '{req_dip}'")
                        print(f"У пациента есть: {snp_marker} = '{pat_gt}' И {star_gene} = '{pat_dip}'")
                        print("----------------------\n")

                    if pat_gt == req_gt and pat_dip == req_dip:
                        reports.append(self._build_report_row(sample, rule, drug, 'COMBINED_STAR_SNP', f"{star_gene}+{snp_marker}", f"{pat_dip}+{pat_gt}"))
                        processed_drugs.add(drug)

            # --- ЭТАП 2: Одиночные правила ---

            # 2.1 Одиночные Star правила
            for _, rule in self.rules_star.iterrows():
                drug = rule['drug']
                gene = rule['gene']
                if drug not in processed_drugs:
                    if gene in star_calls and star_calls[gene] == rule['diplotype']:
                        reports.append(self._build_report_row(sample, rule, drug, 'STAR', gene, star_calls[gene], rule.get('phenotype', '')))
                        processed_drugs.add(drug)

            # 2.2 Одиночные SNP правила
            for _, rule in self.rules_snp.iterrows():
                drug = rule['drug']
                marker = rule['marker_id']
                if drug not in processed_drugs:
                    if marker in snp_calls and snp_calls[marker] == rule['genotype_rule']:
                        reports.append(self._build_report_row(sample, rule, drug, 'SNP', marker, snp_calls[marker]))
                        processed_drugs.add(drug)

        return pd.DataFrame(reports)

    def _build_report_row(self, sample, rule_row, drug, rule_type, target, result, phenotype=''):
        return {
            'sample': sample,
            'drug': drug,
            'rule_type': rule_type,
            'target': target,
            'patient_result': result,
            'phenotype': phenotype,
            'effect': rule_row.get('effect', ''),
            'recommendation': rule_row.get('recommendation', '')
        }
# ==========================================
# 4. REPORTER MODULE (HTML & Summary)
# ==========================================
class Reporter:
    def generate_html(self, report_df, output_path):
        html_content = """
        <!DOCTYPE html>
        <html lang="ru">
        <head>
            <meta charset="UTF-8">
            <title>Фармакогенетический профиль пациента</title>
            <style>
                /* Общие настройки страницы */
                body {
                    font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                    background-color: #f4f7f6; /* Мягкий фон, не режет глаза */
                    color: #333333;
                    margin: 40px auto;
                    max-width: 1200px;
                    padding: 0 20px;
                    line-height: 1.5;
                }
                
                /* Заголовки */
                h1 {
                    color: #2c3e50;
                    text-align: center;
                    border-bottom: 2px solid #3498db;
                    padding-bottom: 15px;
                    margin-bottom: 40px;
                    text-transform: uppercase;
                    letter-spacing: 1px;
                }
                h2 {
                    color: #2c3e50;
                    margin-top: 40px;
                    font-size: 1.2em;
                    background: #ffffff;
                    padding: 12px 20px;
                    border-left: 5px solid #3498db;
                    border-radius: 6px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
                }

                /* Контейнер для таблицы с закругленными краями и тенью */
                .table-container {
                    background: #ffffff;
                    border-radius: 8px;
                    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
                    overflow: hidden; /* Важно для закругления углов у таблицы */
                    margin-bottom: 40px;
                }

                /* Стилизация самой таблицы */
                table {
                    width: 100%;
                    border-collapse: collapse;
                    text-align: left;
                }
                th, td {
                    padding: 16px 20px;
                }
                th {
                    background-color: #2c3e50; /* Темно-синий/серый цвет шапки */
                    color: #ffffff;
                    font-weight: 600;
                    font-size: 0.9em;
                    text-transform: uppercase; /* Все заголовки в верхнем регистре */
                    letter-spacing: 0.05em;
                }
                
                /* Чередование цветов строк для удобства чтения */
                tr {
                    border-bottom: 1px solid #edf2f7;
                }
                tr:last-child {
                    border-bottom: none;
                }
                tr:nth-child(even) {
                    background-color: #f8fafc; /* Очень светлый серый для четных строк */
                }
                tr:hover {
                    background-color: #e2e8f0; /* Подсветка строки при наведении */
                    transition: background-color 0.2s ease;
                }

                /* Декоративные элементы (Бейджи) */
                .badge {
                    display: inline-block;
                    padding: 4px 10px;
                    border-radius: 12px;
                    background-color: #e1effe;
                    color: #1e429f;
                    font-weight: 600;
                    font-size: 0.85em;
                    text-align: center;
                    min-width: 40px;
                }
                
                /* Нормализация регистра */
                .drug-name {
                    text-transform: lowercase; /* Принудительно в нижний регистр */
                    font-weight: bold;
                    color: #2c3e50;
                }
                
                .empty-msg {
                    text-align: center;
                    color: #7f8c8d;
                    font-style: italic;
                    margin-top: 50px;
                }
            </style>
        </head>
        <body>
            <h1>Клинический фармакогенетический отчёт</h1>
        """
        
        if report_df.empty:
            html_content += "<p class='empty-msg'>Нет данных для формирования рекомендаций.</p>"
        else:
            for sample, df in report_df.groupby('sample'):
                html_content += f"<h2>Пациент ID: {sample}</h2>"
                html_content += "<div class='table-container'>"
                html_content += "<table>"
                
                # Шапка таблицы (без колонки confidence)
                html_content += """
                <tr>
                    <th>Препарат</th>
                    <th>Таргет (Ген/Маркер)</th>
                    <th>Результат пациента</th>
                    <th>Фенотип</th>
                    <th>Ожидаемый эффект</th>
                    <th>Клиническая рекомендация</th>
                </tr>
                """
                for _, row in df.iterrows():
                    # Приводим препарат к нижнему регистру на уровне Python + CSS
                    drug_name = str(row['drug']).lower()
                    
                    # Если фенотип пустой, не рисуем синий бейдж
                    phenotype_html = f"<span class='badge'>{row['phenotype']}</span>" if row['phenotype'] else "—"
                    effect = row['effect'] if pd.notna(row['effect']) else "—"
                    
                    html_content += f"""
                    <tr>
                        <td class='drug-name'>{drug_name}</td>
                        <td>{row['target']}</td>
                        <td><b>{row['patient_result']}</b></td>
                        <td>{phenotype_html}</td>
                        <td>{effect}</td>
                        <td>{row['recommendation']}</td>
                    </tr>
                    """
                html_content += "</table></div>"
                
        html_content += "</body></html>"
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
# ==========================================
# MAIN EXECUTION FLOW
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="PGx Pipeline")
    parser.add_argument('--vcf', required=True, help="Path to input VCF")
    parser.add_argument('--rules_dir', required=True, help="Directory with TSV rules")
    parser.add_argument('--out_dir', required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.join(args.out_dir, 'reports'), exist_ok=True)

    print("[1/5] Extracting variants from VCF...")
    extractor = VCFExtractor(os.path.join(args.rules_dir, 'pgx_variants.tsv'))
    extracted_df, samples = extractor.parse_vcf(args.vcf)
    extracted_df.to_csv(os.path.join(args.out_dir, 'extracted_variants.tsv'), sep='\t', index=False)

    print("[2/5] Calling Genotypes and Diplotypes...")
    caller = PGxCaller(os.path.join(args.rules_dir, 'star_alleles.tsv'))
    gene_calls, coverage_stats = caller.call_genes(extracted_df)
    gene_calls.to_csv(os.path.join(args.out_dir, 'gene_calls.tsv'), sep='\t', index=False)
    coverage_stats.to_csv(os.path.join(args.out_dir, 'coverage_report.tsv'), sep='\t', index=False)

    print("[3/5] Applying Clinical Rules...")
    engine = RulesEngine(args.rules_dir)
    drug_report = engine.apply_rules(gene_calls)
    drug_report.to_csv(os.path.join(args.out_dir, 'reports', 'drug_report.tsv'), sep='\t', index=False)

    print("[4/5] Generating Summary...")
    # Summary - агрегация по пациентам
    summary = drug_report.groupby(['sample', 'drug']).size().reset_index(name='rule_matches')
    summary.to_csv(os.path.join(args.out_dir, 'reports', 'summary_by_sample.tsv'), sep='\t', index=False)

    print("[5/5] Generating HTML Report...")
    reporter = Reporter()
    reporter.generate_html(drug_report, os.path.join(args.out_dir, 'reports', 'sample_report.html'))

    print(f"Pipeline finished! Results saved to {args.out_dir}")

if __name__ == "__main__":
    main()
