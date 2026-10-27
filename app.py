# -*- coding: utf-8 -*-
import streamlit as st
import numpy as np
import pandas as pd
import random
import textwrap
import time

st.set_page_config(page_title="Genetik Veri Simülatörü", page_icon="🧬", layout="wide")

def generate_reference_genome(num_chr, chr_len):
    bases = np.array(['A', 'C', 'G', 'T'])
    genome = {}
    fasta_str = ""
    
    for i in range(1, num_chr + 1):
        chr_name = f"Chr{i:02d}"
        sequence_array = np.random.choice(bases, chr_len)
        sequence_str = "".join(sequence_array)
        genome[chr_name] = sequence_str
        
        fasta_str += f">{chr_name} length={chr_len} simulated_for_education\n"
        fasta_str += "\n".join(textwrap.wrap(sequence_str, 80)) + "\n"
        
    return genome, fasta_str

def generate_vcf(genome, num_snps, num_samples, chr_len, gt_probs):
    sample_names = [f"Ogrenci_Ornek_{i+1:02d}" for i in range(num_samples)]
    chromosomes = list(genome.keys())
    
    vcf_chrs = []
    vcf_pos = []
    vcf_refs = []
    vcf_alts = []
    
    for _ in range(num_snps):
        c = random.choice(chromosomes)
        p = random.randint(1, chr_len)
        ref_base = genome[c][p-1]
        alt_base = random.choice([b for b in ['A', 'C', 'G', 'T'] if b != ref_base])
        
        vcf_chrs.append(c)
        vcf_pos.append(p)
        vcf_refs.append(ref_base)
        vcf_alts.append(alt_base)
    
    df_snps = pd.DataFrame({'#CHROM': vcf_chrs, 'POS': vcf_pos, 'REF': vcf_refs, 'ALT': vcf_alts})
    df_snps = df_snps.sort_values(by=['#CHROM', 'POS']).reset_index(drop=True)
    
    df_snps['ID'] = [f"SimSNP_{i}" for i in range(1, num_snps + 1)]
    df_snps['QUAL'] = "99.0"
    df_snps['FILTER'] = "PASS"
    df_snps['INFO'] = "PR"
    df_snps['FORMAT'] = "GT"
    
    df_fixed = df_snps[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
    
    genotypes = ['0/0', '0/1', '1/1', './.']
    gt_matrix = np.random.choice(genotypes, size=(num_snps, num_samples), p=gt_probs)
    
    df_samples = pd.DataFrame(gt_matrix, columns=sample_names)
    df_vcf = pd.concat([df_fixed, df_samples], axis=1)
    
    vcf_header = "##fileformat=VCFv4.2\n"
    vcf_header += "##fileDate=20260902\n"
    vcf_header += "##source=Genetik_Egitim_Simulatoru_Web\n"
    for chr_name in chromosomes:
        vcf_header += f"##contig=<ID={chr_name},length={chr_len}>\n"
    vcf_header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    
    vcf_data = df_vcf.to_csv(sep='\t', index=False)
    return vcf_header + vcf_data, df_vcf

def main():
    st.title("🧬 Eğitim Amaçlı Sentetik Genetik Veri Simülatörü")
    st.markdown("Biyoinformatik dersleri ve uygulamalı eğitimler için **FASTA Referans Genomu** ve buna tam uyumlu **VCF Varyant Dosyası** üreten arayüz.")
    st.divider()

    with st.sidebar:
        st.header("⚙️ Simülasyon Parametreleri")
        
        st.subheader("1. Popülasyon & Genom")
        num_samples = st.number_input("Örnek (Birey) Sayısı", min_value=1, max_value=500, value=50, step=1)
        num_snps = st.number_input("Oluşturulacak SNP Sayısı", min_value=10, max_value=100000, value=5000, step=500)
        num_chr = st.number_input("Kromozom Sayısı", min_value=1, max_value=50, value=7, step=1)
        chr_len = st.number_input("Her Bir Kromozom Uzunluğu (bp)", min_value=1000, max_value=50000000, value=2000000, step=10000)
        
        st.divider()
        st.subheader("2. Genotip Olasılıkları (%)")
        st.markdown("Toplamlarının %100 olması gerekmez, sistem otomatik olarak oranlayacaktır.")
        
        p_ref = st.slider("Homozigot Referans (0/0)", 0, 100, 45)
        p_het = st.slider("Heterozigot (0/1)", 0, 100, 20)
        p_alt = st.slider("Homozigot Varyant (1/1)", 0, 100, 30)
        p_mis = st.slider("Kayıp Veri / Dropout (./.)", 0, 100, 5)
        
        total_p = p_ref + p_het + p_alt + p_mis
        if total_p == 0:
            st.error("Olasılıkların toplamı 0 olamaz!")
            gt_probs = [0.25, 0.25, 0.25, 0.25]
        else:
            gt_probs = [p_ref/total_p, p_het/total_p, p_alt/total_p, p_mis/total_p]
            
        st.info(f"**Güncel Dağılım:**\n\n0/0: %{gt_probs[0]*100:.1f} | 0/1: %{gt_probs[1]*100:.1f}\n\n1/1: %{gt_probs[2]*100:.1f} | ./. : %{gt_probs[3]*100:.1f}")

    col1, col2 = st.columns([1, 1])

    with col1:
        if st.button("🚀 Sentetik Verisetini Oluştur", use_container_width=True, type="primary"):
            with st.spinner("1. Referans Genom (FASTA) üretiliyor..."):
                start_time = time.time()
                genome_dict, fasta_string = generate_reference_genome(num_chr, chr_len)
                
            with st.spinner("2. Genotipler ve VCF matrisi hesaplanıyor..."):
                vcf_string, df_preview = generate_vcf(genome_dict, num_snps, num_samples, chr_len, gt_probs)
                end_time = time.time()
                
            st.session_state['fasta_data'] = fasta_string
            st.session_state['vcf_data'] = vcf_string
            st.session_state['df_preview'] = df_preview
            st.session_state['time_taken'] = round(end_time - start_time, 2)
            st.success(f"✅ İşlem tamamlandı! ({st.session_state['time_taken']} saniye)")

    if 'vcf_data' in st.session_state:
        st.subheader("👀 VCF Önizleme (İlk 15 SNP)")
        st.dataframe(st.session_state['df_preview'].head(15), use_container_width=True)
        
        st.divider()
        st.subheader("📥 Dosyaları Bilgisayara İndir")
        
        dl_col1, dl_col2 = st.columns(2)
        with dl_col1:
            st.download_button(
                label="📁 FASTA Referans Genomunu İndir (.fasta)",
                data=st.session_state['fasta_data'],
                file_name="Egitim_Referans_Genom.fasta",
                mime="text/plain",
                use_container_width=True
            )
        with dl_col2:
            st.download_button(
                label="📊 VCF Genotip Dosyasını İndir (.vcf)",
                data=st.session_state['vcf_data'],
                file_name="Egitim_Sentetik_GBS.vcf",
                mime="text/plain",
                use_container_width=True
            )

if __name__ == "__main__":
    main()
