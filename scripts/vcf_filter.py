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
try:
    vcf_reader = vcf.Reader(open(vcf_file, "r"))  # Ouvre le fichier en mode lecture
except Exception as e:
    print(f"Erreur lors de la lecture du fichier VCF : {e}")
    exit()  # Arrête le script si une erreur survient à la lecture

## ---------Ajout de filtres pour garder uniquement les variants de bonne qualité-------------##

# Définition des seuils fixes selon le protocole
QUAL_THRESHOLD = 30  # Score de qualité minimal requis
AF_THRESHOLD = 0.01  # Fréquence allélique minimale acceptée
DP_THRESHOLD = 10  # Profondeur de lecture minimale

# Initialisation d'une liste pour stocker les variants extraits après filtrage
variants = []

# Parcours de chaque ligne du fichier VCF
for record in vcf_reader:
    # Vérification du score QUAL
    if record.QUAL < QUAL_THRESHOLD:
        continue  # Ignore la mutation si la qualité est trop faible
    
    # Vérification de la fréquence allélique (AF)
    if record.INFO.get("AF", [None])[0] is not None and record.INFO["AF"][0] < AF_THRESHOLD:
        continue  # Ignore la mutation si la fréquence allélique est inférieure au seuil
    
    # Vérification de la profondeur de lecture (DP)
    if "DP" in record.INFO and record.INFO["DP"] < DP_THRESHOLD:
        continue  # Ignore la mutation si la profondeur est insuffisant
    
    # Extraction des informations importantes sur la mutation
    variant_info = {
        "CHROM": record.CHROM,  # Chromosome sur lequel la mutation est située
        "POS": record.POS,  # Position exacte de la mutation
        "REF": record.REF,  # Allèle de référence
        "ALT": ",".join(str(alt) for alt in record.ALT),  # Allèles alternatifs trouvés
        "QUAL": record.QUAL,  # Score de qualité de la mutation
        "FILTER": ",".join(record.FILTER) if record.FILTER else "PASS",  # Indication des filtres appliqués
        "DP": record.INFO.get("DP", "NA"),  # Profondeur de lecture
    }
    
    # Ajout du variant à la liste après filtrage
    variants.append(variant_info)

# Conversion de la liste en tableau `pandas` pour une manipulation plus simple
if variants:
    df = pd.DataFrame(variants)
    # Affichage des statistiques des données extraites
    print(df.describe())  # Donne un aperçu des valeurs numériques des mutations filtrées
    print("Aperçu des variants extraits :")
    print(df.head())  # Affiche les 5 premières lignes des mutations conservées
    print("\nNombre total de variants extraits :", len(df))  # Nombre total de mutations retenues
else:
    print("Aucun variant retenu après filtrage.")
    exit()  # Arrête le script si aucun variant n'est retenu

## --------------------- Sauvegarde dans dossier puis un fichier TSV---------------- ##

# Définition du dossier de sauvegarde pour stocker les résultats filtrés
base_dir = os.path.dirname(os.path.abspath(vcf_file))  # Récupérer le dossier contenant le fichier VCF
output_dir = os.path.join(base_dir, "filter_variant")  # Créer le sous-dossier filter_variant dans data
os.makedirs(output_dir, exist_ok=True)  # Crée `data/filter_variant/` s'il n'existe pas

# Sauvegarde des résultats filtrés dans un fichier TSV
output_file = os.path.join(output_dir, os.path.splitext(os.path.basename(vcf_file))[0] + "_filtered.tsv")
df.to_csv(output_file, index=False, sep="\t")  # Export du fichier TSV avec tabulation

# Affichage du chemin du fichier de sortie
print(f"\n Résultats sauvegardés dans : {output_file}")
