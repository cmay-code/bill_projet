#!/bin/bash

# Vérification du dossier donné en argument
if [ -z "$1" ]; then
    echo " Erreur : Aucun dossier de données spécifié !"
    echo " Utilisation : bash fusion_vcf_snp.sh /chemin/vers/dossier_vcf"
    exit 1
fi

data_dir=$1  # Dossier contenant les fichiers VCF
fusion_dir="${data_dir}/fusion_snp"  # Dossier où stocker les fichiers fusionnés
mkdir -p "$fusion_dir"  # Création du dossier fusion s'il n'existe pas

# Fusion pour chaque passage
for passage in P30 P65 P90; do
    echo " Fusion des fichiers SNP pour $passage..."

    # Fusion des fichiers froids (1 à 5)
    first_file=$(ls ${data_dir}/${passage}/${passage}-1.trimed1000.snp.vcf | head -n 1)
    if [ -f "$first_file" ]; then
        grep "^#" "$first_file" > ${fusion_dir}/${passage}_cold_snp.vcf  # Garder l'entête
        cat ${data_dir}/${passage}/${passage}-{1..5}.trimed1000.snp.vcf | grep -v "^#" >> ${fusion_dir}/${passage}_cold_snp.vcf
        echo " Fusion terminée pour ${passage}_cold_snp.vcf"
    else
        echo "⚠️ Aucun fichier froid trouvé pour ${passage}, fusion ignorée."
    fi

    # Fusion des fichiers chauds (6 à 10)
    first_file=$(ls ${data_dir}/${passage}/${passage}-6.trimed1000.snp.vcf | head -n 1)
    if [ -f "$first_file" ]; then
        grep "^#" "$first_file" > ${fusion_dir}/${passage}_hot_snp.vcf  # Garder l'entête
        cat ${data_dir}/${passage}/${passage}-{6..10}.trimed1000.snp.vcf | grep -v "^#" >> ${fusion_dir}/${passage}_hot_snp.vcf
        echo " Fusion terminée pour ${passage}_hot_snp.vcf"
    else
        echo " Aucun fichier chaud trouvé pour ${passage}, fusion ignorée."
    fi
done

echo " Toutes les fusions des SNP sont terminées ! Les fichiers sont dans ${fusion_dir}/"
