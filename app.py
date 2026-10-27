import streamlit as st
import pandas as pd
import numpy as np
import io

st.set_page_config(page_title="VCF Veri Üretici", layout="centered")

st.title("🧬 Sentetik VCF-Excel Veri Üretici")
st.markdown("Bu program, test süreçlerinde kullanılmak üzere istediğiniz boyutlarda yapay bir VCF-Excel dosyası oluşturur.")

# Kullanıcı Girişleri
col1, col2 = st.columns(2)
with col1:
    n_genotypes = st.number_input("Toplam Bitki (Genotip) Sayısı:", min_value=1, max_value=5000, value=100)
with col2:
    n_snps = st.number_input("Toplam SNP Sayısı:", min_value=1, max_value=50000, value=1000)

st.sidebar.header("Genetik Ayarlar")
het_rate = st.sidebar.slider("Heterozigotluk Oranı (0/1):", 0.0, 0.5, 0.1)
alt_rate = st.sidebar.slider("Mutasyon Oranı (1/1):", 0.0, 0.5, 0.2)
ref_rate = 1.0 - (het_rate + alt_rate)

if st.button("VCF-Excel Dosyasını Oluştur"):
    with st.spinner("Veri seti oluşturuluyor..."):
        # 1. Genotip ve SNP isimlerini hazırla
        samples = [f"Variety_{i+1}" for i in range(int(n_genotypes))]
        snp_names = [f"SNP_{i+1}" for i in range(int(n_snps))]
        
        # 2. VCF Bilgi Sütunlarını Oluştur
        df_info = pd.DataFrame({
            'CHROM': ['chr1'] * int(n_snps),
            'POS': [i * 100 for i in range(1, int(n_snps) + 1)],
            'ID': snp_names,
            'REF': ['A'] * int(n_snps),
            'ALT': ['G'] * int(n_snps),
            'QUAL': [100] * int(n_snps),
            'FILTER': ['PASS'] * int(n_snps),
            'INFO': ['.'] * int(n_snps),
            'FORMAT': ['GT'] * int(n_snps)
        })
        
        # 3. Rastgele Genotip Verilerini Üret
        calls = np.random.choice(["0/0", "0/1", "1/1"], 
                                 size=(int(n_snps), int(n_genotypes)), 
                                 p=[ref_rate, het_rate, alt_rate])
        
        df_genotypes = pd.DataFrame(calls, columns=samples)
        
        # 4. Birleştir
        full_vcf = pd.concat([df_info, df_genotypes], axis=1)
        
        # 5. Excel'e Dönüştür
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            full_vcf.to_excel(writer, index=False, sheet_name='VCF_Data')
        
        st.success(f"✅ {n_genotypes} bitki ve {n_snps} SNP içeren dosya hazır!")
        
        st.download_button(
            label="📥 Oluşturulan VCF-Excel Dosyasını İndir",
            data=output.getvalue(),
            file_name=f"sentetik_vcf_{n_genotypes}v_{n_snps}s.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )