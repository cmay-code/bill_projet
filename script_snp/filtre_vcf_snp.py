import os
import sys
import pandas as pd

# Dossiers d'entrée et de sortie
fusion_dir = sys.argv[1]  # Dossier contenant les fichiers fusionnés
output_dir = os.path.join(fusion_dir, "filtre")  # Dossier où stocker les fichiers filtrés
os.makedirs(output_dir, exist_ok=True)  # Crée le dossier s'il n'existe pas

# Filtrage des fichiers fusionnés
for file in os.listdir(fusion_dir):
    if file.endswith(".vcf"):  # On ne traite que les fichiers VCF
        vcf_path = os.path.join(fusion_dir, file)
        output_csv = os.path.join(output_dir, file.replace(".vcf", "_filtre.csv"))

        # Lire le VCF et extraire les données utiles
        filtered_variants = []
        with open(vcf_path, "r") as vcf:
            for line in vcf:
                if line.startswith("#"):  # Ignorer les commentaires
                    continue

                cols = line.strip().split("\t")
                _, pos, id_, ref, alt, qual, filt, info, *others = cols  # Ignore CHROM et autres colonnes inutiles

                # Extraction des valeurs de INFO
                info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in info.split(";") if "=" in kv}

                # Récupérer les valeurs nécessaires avec valeurs par défaut
                af = float(info_dict.get("AF", 0))  # Fréquence allélique
                dp = int(info_dict.get("DP", 0))  # Profondeur de lecture

                # Appliquer les filtres réadaptés pour SNPs
                if float(qual) > 30 and filt == "PASS" and af >= 0.1 and dp >= 10:
                    filtered_variants.append([pos, id_, ref, alt, qual, filt, info])  # Supprime les FIELD_X

        # Récupérer les noms de colonnes essentielles
        header = ["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

        # Convertir en CSV
        df = pd.DataFrame(filtered_variants, columns=header)
        df.to_csv(output_csv, sep="\t", index=False)

        print(f"✅ Fichier filtré enregistré : {output_csv}")

print(" Tous les fichiers VCF fusionnés ont été filtrés et enregistrés dans le dossier `filtre/` !")
