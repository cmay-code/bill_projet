#!/bin/bash

# ------------------------------------------ #
# SCRIPT DE FUSION DES FICHIERS SNP & SV    #
# ------------------------------------------ #
# Ce script fusionne les fichiers SNP et SV 
# pour chaque passage (P15, P50, P90) en veillant 
# à préserver l’en-tête des fichiers VCF.
# ------------------------------------------ #

#  Définition du dossier contenant les données (chemin automatique)
data_dir="$(pwd)/../data"

# Vérification de l'existence du dossier contenant les fichiers à fusionner
if [ ! -d "$data_dir" ]; then
    echo " Erreur : Le dossier 'data' n'existe pas dans $(pwd). Vérifie ton chemin."
    exit 1
fi

#  Liste des passages expérimentaux à traiter
for passage in P15 P50 P90; do
    echo " Fusion des fichiers SNP et SV pour ${passage}..."
    dossier="${data_dir}/${passage}"

    # Vérifier si le dossier du passage existe
    if [ ! -d "${dossier}" ]; then
        echo " Erreur : Le dossier ${dossier} n'existe pas !"
        continue
    fi

    # Fusionner les fichiers SNP et SV pour chaque index (1 à 10)
    for i in {1..10}; do
        snp_file="${dossier}/${passage}-${i}.trimed1000.snp.vcf"
        sv_file="${dossier}/${passage}-${i}.trimed1000.sv_sniffles.vcf"
        output_file="${dossier}/${passage}_${i}.vcf"

        # Vérification et fusion des fichiers SNP et SV
        if [[ -f "$snp_file" && -f "$sv_file" ]]; then
            echo " Fusion en cours pour ${output_file}..."
            grep "^#" "$snp_file" > "$output_file"  # Récupérer l’en-tête du SNP
            (grep -v "^#" "$snp_file"; grep -v "^#" "$sv_file") | sort -u >> "$output_file"  # Fusionner sans l'en-tête
            echo " Fusionnée avec succès : $output_file"

        elif [[ -f "$snp_file" ]]; then
            cp "$snp_file" "$output_file"
            echo "Copie seule du SNP (pas de SV trouvé) : $snp_file → $output_file"

        elif [[ -f "$sv_file" ]]; then
            cp "$sv_file" "$output_file"
            echo " Copie seule du SV (pas de SNP trouvé) : $sv_file → $output_file"

        else
            echo " Aucun fichier SNP/SV trouvé pour ${passage} index ${i} !"
        fi
    done
done

echo " Fusion SNP + SV terminée avec succès pour tous les passages !"
