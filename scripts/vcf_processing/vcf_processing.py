## ----------------------Importation des librairies nécessaires-------------##

import vcf  # Permet la lecture des fichiers VCF
import pandas as pd  # Permet la manipulation des données sous forme de tableau
import argparse  # Pour récupérer les arguments en ligne de commande
import os  # Gestion des fichiers et des chemins


## -----------------------Lecture et extraction des données------------------##

# Définition du parser pour récupérer le fichier en argument
parser = argparse.ArgumentParser(description="Analyse d'un fichier VCF")
parser.add_argument("vcf_file", help="Chemin du fichier VCF à analyser")  # Argument obligatoire : fichier VCF
args = parser.parse_args()  # Lecture des arguments fournis par l'utilisateur

# Récupération du chemin du fichier VCF entré par l'utilisateur
vcf_file = args.vcf_file

# Vérification que le fichier existe avant de l'ouvrir
if not os.path.exists(vcf_file):  
    print(f"Erreur : Le fichier '{vcf_file}' n'existe pas. Vérifiez le chemin.")
    exit()  # Arrête le script si le fichier est introuvable

# Tentative d'ouverture et de lecture du fichier VCF
try: #utilisé pour gérer les erreurs lorsqu'on ouvre le fichier VCF.
    vcf_reader = vcf.Reader(open(vcf_file, "r"))  # Ouvre le fichier en mode lecture
except Exception as e: # Capture toutes les erreurs possibles lors de l'ouverture du fichier VCF.
# Si une erreur survient (fichier introuvable, corrompu, mauvais format, etc.),
# elle est stockée dans la variable "e" et un message d'erreur sera affiché.
    print(f"Erreur lors de la lecture du fichier VCF : {e}")
    exit()  # Arrête le script si une erreur survient à la lecture



## ---------Ajoute de filtres pour garder uniquement les variants de bonne qualité-------------##



# Initialisation d'une liste pour stocker les variants extraits
variants = []

# Parcours de chaque ligne du fichier VCF
for record in vcf_reader:
    # Extraction des informations importantes
    variant_info = {
        "CHROM": record.CHROM,  # Chromosome
        "POS": record.POS,  # Position du variant sur le chromosome
        "REF": record.REF,  # Allèle de référence
        "ALT": ",".join(str(alt) for alt in record.ALT),  # Allèles alternatifs (peut y en avoir plusieurs), concaténation en chaîne séparé par ,
        "QUAL": record.QUAL,  # Score de qualité du variant
        "FILTER": ",".join(record.FILTER) if record.FILTER else "PASS",  # Filtre appliqué au variant
    }

    # Application des filtres
    if record.INFO.get("AF", [None])[0] is not None and (record.INFO["AF"][0] < 0.01 or record.INFO["AF"][0] > 0.99):
        continue  # On ignore ce variant si sa fréquence allélique est trop basse/élevée

    if "DP" in record.INFO and record.INFO["DP"] < 10:
        continue  # On ignore ce variant si la profondeur de lecture est faible

    if "SUPPORT" in record.INFO and record.INFO["SUPPORT"] < 10:
        continue  # On ignore ce variant si le support est faible (SV)

    if "SVLEN" in record.INFO and record.INFO["SVLEN"] < 50:
        continue  # On ignore ce variant si la taille est inférieure à 50 bp (SV)

    if record.FILTER and "PASS" not in record.FILTER:
        continue  # On ignore ce variant si le filtre n'est pas "PASS"

    # Ajout du variant à la liste
    variants.append(variant_info)

# Conversion de la liste en tableau `pandas` pour une manipulation plus simple
df = pd.DataFrame(variants)

#Affichage d'un résumé claire.
print(df.describe())  # Donne un aperçu des valeurs numériques

# Affichage des 5 premières lignes du tableau pour vérifier les données extraites
print("Aperçu des variants extraits :")
print(df.head())

# Affichage d'un résumé des données
print("\n Nombre total de variants extraits :", len(df))
print("\n Colonnes disponibles dans le tableau :", df.columns)




## --------------------- Sauvegarde dans dossier puis un fichier TSV---------------- ##

# Définition du chemin du dossier "bill_projet"
project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Remonte au dossier bill_projet

# Définition du dossier de sauvegarde "new_variants" dans "bill_projet"
output_dir = os.path.join(project_dir, "new_variants")

# Vérification et création du dossier s'il n'existe pas
if not os.path.exists(output_dir):
    os.makedirs(output_dir)  # Crée le dossier si nécessaire

# Sauvegarde des résultats filtrés dans le bon dossier
output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(vcf_file))[0] + "_filtered.tsv")
df.to_csv(output_file, index=False, sep="\t")  # Export avec tabulation

print(f"\n✅ Résultats sauvegardés dans : {output_file}")

