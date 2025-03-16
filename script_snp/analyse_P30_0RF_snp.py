import os
import sys
import pandas as pd
from Bio import SeqIO
import re

print(" Initialisation du script...")

# 1️⃣ Définition du dossier contenant les fichiers filtrés
if len(sys.argv) > 1:
    dossier_filtrage = sys.argv[1]  # Chemin passé en argument
else:
    dossier_filtrage = "./filtrage"  # Valeur par défaut

# Vérifier si le dossier existe
if not os.path.exists(dossier_filtrage):
    print(f" Erreur : Le dossier {dossier_filtrage} n'existe pas !")
    exit()

print(f" Lecture des fichiers dans le dossier : {dossier_filtrage}")

# 2️⃣ Lecture des fichiers de mutations filtrées
fichiers = [f for f in os.listdir(dossier_filtrage) if f.endswith(".csv")]
if not fichiers:
    print(" Aucun fichier CSV trouvé dans le dossier filtrage !")
    exit()

# Définir le dossier de sortie dans le dossier de filtrage
dossier_sortie = os.path.join(dossier_filtrage, "mutation_P30_Orf")

# ✅ Vérifier et créer le dossier de sortie
if not os.path.exists(dossier_sortie):
    os.makedirs(dossier_sortie)
    print(f"📁 Dossier créé : {dossier_sortie}")

# Stockage des mutations
mutations = {"P30": {"hot": [], "cold": []}}

for fichier in fichiers:
    chemin_fichier = os.path.join(dossier_filtrage, fichier)
    print(f" Lecture de {fichier}...")

    try:
        df = pd.read_csv(chemin_fichier, sep="\t")
        print(f"Colonnes disponibles dans {fichier} : {df.columns.tolist()}")
    except Exception as e:
        print(f" Erreur de lecture pour {fichier} : {e}")
        continue

    # Vérification de la colonne INFO
    if "INFO" not in df.columns:
        print(f"Fichier {fichier} ignoré : colonne INFO manquante !")
        continue

    # Détermination de la condition hot/cold
    condition = "hot" if "hot" in fichier else "cold"

    # Stocker les mutations de P30
    for _, row in df.iterrows():
        try:
            info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in row["INFO"].split(";") if "=" in kv}

            variant = {
                "pos": int(row["POS"]),
                "alt": row["ALT"],
                "ref": row["REF"],
                "coverage": int(info_dict.get("DP", "0")),
                "Af": float(info_dict.get("AF", "0.0"))
            }
            mutations["P30"][condition].append(variant)

        except Exception as e:
            print(f"⚠️ Erreur lors de l'extraction des données pour une ligne dans {fichier} : {e}")
            continue    

print("✅ Lecture des fichiers terminée !")

# 3️⃣ Lecture du fichier ORF.fasta et annotation des mutations
orf_data = []
fichier_orf = input("Entrez le chemin du fichier ORF.fasta : ")

if not os.path.exists(fichier_orf):
    print("❌ Erreur : Le fichier ORF.fasta n'existe pas.")
    exit()

for record in SeqIO.parse(fichier_orf, "fasta"):
    desc = record.description
    match = re.search(r'location=(?:complement\()?([0-9]+)\.\.([0-9]+)', desc)
    if match:
        orf_start = int(match.group(1))
        orf_end = int(match.group(2))
        orf_name = re.search(r'\[protein=(.+?)\]', desc)
        orf_name = orf_name.group(1) if orf_name else "Unknown_ORF"
        orf_seq = str(record.seq)
        orf_data.append({"name": orf_name, "start": orf_start, "end": orf_end, "sequence": orf_seq})
    else:
        print(f"⚠️ Erreur : Position non trouvée pour {desc}")

print("✅ Lecture des ORF terminée !")

# 4️⃣ Vérification des mutations dans les ORF et récupération de la base mutée
def check_orf(mutations, orf_data):
    annotated_mutations = []
    for mutation in mutations:
        mutation["ORF"] = "None"
        mutation["ORF_start"] = None
        mutation["ORF_end"] = None
        mutation["Base_ORF"] = "None"

        mutation_pos = mutation["pos"]
        found = False

        for orf in orf_data:
            orf_start, orf_end = orf["start"], orf["end"]

            if orf_start <= mutation_pos <= orf_end:
                mutation["ORF"] = orf["name"]
                mutation["ORF_start"] = orf_start
                mutation["ORF_end"] = orf_end
                mutation["Base_ORF"] = orf["sequence"][mutation_pos - orf_start]
                found = True
                break
        
        if found:
            annotated_mutations.append(mutation)
        else:
            print(f"⚠️ Mutation à la position {mutation_pos} n'a pas trouvé d'ORF correspondante.")
    
    return annotated_mutations

mutations["P30"]["hot"] = check_orf(mutations["P30"]["hot"], orf_data)
mutations["P30"]["cold"] = check_orf(mutations["P30"]["cold"], orf_data)

# 5️⃣ Sauvegarde des fichiers finaux
df_hot = pd.DataFrame(mutations["P30"]["hot"])
df_cold = pd.DataFrame(mutations["P30"]["cold"])

df_hot.to_csv(os.path.join(dossier_sortie, "mutations_annotées_P30_hot.csv"), index=False)
df_cold.to_csv(os.path.join(dossier_sortie, "mutations_annotées_P30_cold.csv"), index=False)

print("✅ Analyse terminée et fichiers enregistrés !")
