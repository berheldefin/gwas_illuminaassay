#######extract_rsID_from_illumina_one
import re
import csv

def clean_rsID(rsID):
    return rsID.strip()

input_file = '/home/berheldefin/new/GSA-24v3-0_A2.csv'
output_file = '/home/berheldefin/new/final/rsID1.csv'

# rsID'leri depolamak için bir küme (set) oluşturuyoruz
rsID_set = set()

# CSV dosyasını okuyup rsID'leri küme içine ekliyoruz (tekrar edenler sadece bir kez eklenir)
with open(input_file, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        for item in row:
            rsID_matches = re.findall(r'rs[0-9]+', item)
            for rsID in rsID_matches:
                rsID_cleaned = clean_rsID(rsID)
                rsID_set.add(rsID_cleaned)

# rsID'leri yeni bir CSV dosyasına yazıyoruz
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['rsID'])  # Başlık olarak "rsID" ekliyoruz
    for rsID in rsID_set:
        writer.writerow([rsID])

print(f"{len(rsID_set)} farklı rsID başarıyla çıktı dosyasına yazıldı.")

##########annotation_two
import pandas as pd

def main():
    # A.tsv dosyasını yükleme
    df_a = pd.read_csv('/home/berheldefin/new/gwas_catalog_v1.0-associations_e110_r2023-07-20.tsv', sep='\t', low_memory=False)

    # B.csv dosyasını yükleme ve sadece rsID sütununu seçme
    df_b = pd.read_csv('/home/berheldefin/new/final/rsID1.csv')

    # Ortak rsID'leri bulma
    common_rsids = df_a[df_a['SNPS'].isin(df_b['rsID'])]

    # Ortak rsID'lerle birleştirme
    merged_df = pd.merge(common_rsids, df_b, left_on='SNPS', right_on='rsID', how='inner')

    # Yeni bir csv dosyasına yazma
    merged_df.to_csv('/home/berheldefin/new/final/annotated.csv', sep="\t", index=False)

    print("Ortak veriler başarıyla annotated.CSV dosyasına kaydedildi.")

if __name__ == '__main__':
    main()

########venn_three
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib_venn import venn2

def main():
    # A.tsv dosyasını yükleme
    df_a = pd.read_csv('/home/berheldefin/new/gwas_catalog_v1.0-associations_e110_r2023-07-20.tsv', sep='\t', low_memory=False)

    # B.csv dosyasını yükleme ve sadece rsID sütununu seçme
    df_b = pd.read_csv('/home/berheldefin/new/final/rsID1.csv')['rsID']

    # Ortak rsID'leri bulma
    common_rsids = df_a[df_a['SNPS'].isin(df_b)]

    # Tekrar eden rsID'leri bir kez sayma
    rsid_counts = Counter(common_rsids['SNPS'])

    # Tekrar eden rsID'leri tek bir kere katarak yeniden DataFrame oluşturma
    unique_common_rsids = pd.DataFrame(list(rsid_counts.items()), columns=['SNPS', 'Count'])

    # Venn diyagramı için verileri hazırlama
    set_a = set(df_a['SNPS'])
    set_b = set(df_b)
    set_common = set(unique_common_rsids['SNPS'])

    # Venn diyagramını çizdirme
    venn2(subsets=(len(set_a - set_common), len(set_b - set_common), len(set_common)),
          set_labels=('gwas catalog', 'illumina product files'))

    # Venn diyagramını gösterme
    plt.title("Venn Diagram of gwas and illumina rsIDs")
    plt.savefig('/home/berheldefin/new/final/venn_diagram.png', format='png')
    plt.show()

    print("Venn_diagram.png başarıyla kaydedildi.")

if __name__ == '__main__':
    main()

########confirmation_four
# gwas rsIDs count
import pandas as pd
from collections import Counter

def main():
    # TSV dosyasını yükleme
    file_path = '/home/berheldefin/new/gwas_catalog_v1.0-associations_e110_r2023-07-20.tsv'
    df = pd.read_csv(file_path, sep='\t', low_memory=False)

    # "SNPS" sütunundaki rsID'leri sayma
    rsid_counts = Counter(df['SNPS'])

    # Tekrar eden rsID'leri sadece bir kez sayacak DataFrame oluşturma
    unique_rsid_counts = pd.DataFrame(rsid_counts.items(), columns=['rsID', 'Count'])

    # Sonuçları yazdırma
    print(f"{len(unique_rsid_counts)} kadar unique gwas rsID bulundu")

if __name__ == '__main__':
    main()


# gwas and illumina common rsIDs count

import pandas as pd

def main():
    # A.tsv dosyasını yükleme
    df_a = pd.read_csv('/home/berheldefin/new/gwas_catalog_v1.0-associations_e110_r2023-07-20.tsv', sep='\t', low_memory=False)

    # B.csv dosyasını yükleme ve sadece rsID sütununu seçme
    df_b = pd.read_csv('/home/berheldefin/new/final/rsID1.csv')

    # Ortak rsID'leri bulma
    common_rsids = df_a[df_a['SNPS'].isin(df_b['rsID'])]

    # Ortak rsID'lerle birleştirme
    merged_df = pd.merge(common_rsids, df_b, left_on='SNPS', right_on='rsID', how='inner')

    # Tekrar eden rsID'leri sadece bir kez tutma
    merged_df = merged_df.drop_duplicates(subset='SNPS', keep='first')

    # Yeni bir csv dosyasına yazma
    merged_df.to_csv('/home/berheldefin/new/final/annotatedunique.csv', sep="\t", index=False)

    print("Ortak veriler başarıyla annotatedunique.csv dosyasına kaydedildi.")

if __name__ == '__main__':
    main()