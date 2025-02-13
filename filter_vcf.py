import os
import sys
import pandas as pd

#  Dossiers d'entrée et de sortie
fusion_dir = sys.argv[1]   # Dossier contenant les fichiers fusionnés
output_dir = os.path.join(fusion_dir, "filtre")  # Dossier où stocker les fichiers filtrés
os.makedirs(output_dir, exist_ok=True)  # Crée le dossier s'il n'existe pas

#  Filtrage des fichiers fusionnés
for file in os.listdir(fusion_dir):
    if file.endswith(".vcf"):  # On ne traite que les fichiers VCF
        vcf_path = os.path.join(fusion_dir, file)
        output_csv = os.path.join(output_dir, file.replace(".vcf", "_filtre.csv"))

        #  Lire le VCF et extraire les données utiles
        filtered_variants = []
        with open(vcf_path, "r") as vcf:
            for line in vcf:
                if line.startswith("#"):  # Ignorer les commentaires
                    continue

                cols = line.strip().split("\t")
                chrom, pos, id_, ref, alt, qual, filt, info = cols[:8]

                #  Extraction des valeurs de INFO
                info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in info.split(";") if "=" in kv}

                # Récupérer les valeurs nécessaires avec valeurs par défaut
                svlen = int(info_dict.get("SVLEN", 0))
                af = float(info_dict.get("AF", 0))
                support = int(info_dict.get("SUPPORT", 0))

                #  Appliquer les filtres
                if float(qual) > 30 and filt == "PASS" and af >= 0.1 and svlen > 10 and support > 10:
                    filtered_variants.append(cols[1:])  # On supprime seulement la colonne CHROM

        #  Convertir en CSV en gardant toutes les colonnes sauf CHROM
        df = pd.DataFrame(filtered_variants, columns=["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
        df.to_csv(output_csv, sep="\t", index=False)

        print(f" Fichier filtré enregistré : {output_csv}")

print(" Tous les fichiers VCF fusionnés ont été filtrés et convertis en CSV dans le dossier `filtre/` !")
