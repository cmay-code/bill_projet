#!/bin/bash

# ------------------------------------------ #
# SCRIPT DE FUSION DES FICHIERS VCF PAR CONDITION (FROID/CHAUD)
# ------------------------------------------ #
# Ce script fusionne les fichiers VCF pour chaque passage expérimental (P15, P50, P90)
# en fonction des conditions expérimentales : Froid (1 à 5) et Chaud (6 à 10).
# Il utilise vcf-concat (VCFtools) pour fusionner les fichiers.
# ------------------------------------------ #

# Vérification de l'argument (chemin du dossier de données)
if [ -z "$1" ]; then
    echo " Erreur : Aucun dossier de données spécifié !"
    echo " Utilisation : bash vcf_fusion.sh /chemin/vers/data"
    exit 1
fi

data_dir=$1  # Le premier argument donné au script correspond au dossier contenant les fichiers VCF

# Vérification de l'existence du dossier contenant les fichiers VCF
if [ ! -d "$data_dir" ]; then
    echo " Erreur : Le dossier spécifié n'existe pas !"
    exit 1
fi

# Vérification que vcftools et vcf-concat sont bien installés
if ! command -v vcf-concat &> /dev/null; then
    echo " Erreur : vcf-concat n'est pas installé."
    echo " Installez-le avec : sudo apt install vcftools"
    exit 1
fi

# Création du dossier de sortie pour les fichiers fusionnés
output_dir="${data_dir}/fusion_results"
mkdir -p $output_dir  # Crée le dossier s'il n'existe pas

# Liste des passages à traiter
for passage in P15 P50 P90; do
    echo ""
    echo " Traitement du passage ${passage}..."
    
    dossier="${data_dir}/${passage}"  # Définition du chemin du passage

    # Vérifier si le dossier du passage existe avant de continuer
    if [ ! -d "$dossier" ]; then
        echo " Erreur : Le dossier ${dossier} n'existe pas !"
        continue
    fi

    # ------------------------------------------ #
    # Fusion des fichiers "Froid" (de 1 à 5)
    # ------------------------------------------ #
    echo " Fusion des fichiers pour ${passage} - Froid..."
    
    vcf-concat ${dossier}/${passage}_1.vcf \
               ${dossier}/${passage}_2.vcf \
               ${dossier}/${passage}_3.vcf \
               ${dossier}/${passage}_4.vcf \
               ${dossier}/${passage}_5.vcf > ${output_dir}/${passage}_froid.vcf

    echo " Fusion terminée pour ${passage} - Froid : ${output_dir}/${passage}_froid.vcf"

    # ------------------------------------------ #
    # Fusion des fichiers "Chaud" (de 6 à 10)
    # ------------------------------------------ #
    echo " Fusion des fichiers pour ${passage} - Chaud..."
    
    vcf-concat ${dossier}/${passage}_6.vcf \
               ${dossier}/${passage}_7.vcf \
               ${dossier}/${passage}_8.vcf \
               ${dossier}/${passage}_9.vcf \
               ${dossier}/${passage}_10.vcf > ${output_dir}/${passage}_chaud.vcf

    echo " Fusion terminée pour ${passage} - Chaud : ${output_dir}/${passage}_chaud.vcf"
    
done

#  Fin du script
echo "--------------------------------------"
echo " Fusion terminée avec succès pour tous les passages !"
