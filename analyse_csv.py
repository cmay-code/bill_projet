import os
import sys
import pandas as pd
from Bio import SeqIO
import re


print(" Initialisation du script...")

# 1️⃣ Définition du dossier où sont stockés les fichiers filtrés
if len(sys.argv) > 1:
    dossier_filtrage = sys.argv[1]  # Récupérer le chemin du dossier en argument
else:
    dossier_filtrage = "./filtrage"  # Valeur par défaut si aucun argument n'est donné

# Vérifier si le dossier existe
if not os.path.exists(dossier_filtrage):
    print(f"Erreur : Le dossier {dossier_filtrage} n'existe pas !")
    exit()

print(f" Lecture des fichiers dans le dossier : {dossier_filtrage}")


# 2️⃣ Lecture des fichiers
fichiers = [f for f in os.listdir(dossier_filtrage) if f.endswith(".csv")] #os.listdir(dossier_filtrage) permet de lister tous les fichiers présents dans le dossier spécifié (dossier_filtrage).
if not fichiers:
    print(" Aucun fichier CSV trouvé dans le dossier filtrage !")
    exit()

mutations = {}
conditions = ["hot", "cold"]

# Lire et stocker les mutations
#Dans cette boucle, on lit les fichier et on les stockes dans le dico mutation de cette manière : 
#mutations = {
    #"P50": {
    #    "hot": [list de mutations à P50 sous hot],
    #   "cold": [list de mutations à P50 sous cold]
    #},
    #"P90": {
    #    "hot": [list de mutations à P90 sous hot],
    #   "cold": [list de mutations à P90 sous cold]
    #}
#}

for fichier in fichiers:
    chemin_fichier = os.path.join(dossier_filtrage, fichier)
    print(f" Lecture de {fichier}...")
    #après la lecture en boucle du fichier, on utilise pandas stocker en Dataframe (df)
    #expect (fonction de panda) qui perment de faire une gestion des erreurs de lecture
    #si ce fichier ne peut pas être lu à cause d'un soucis, le script ne s'arrête pas mais
    #il passe au fichier suivant grâce à continue
    try:
        df = pd.read_csv(chemin_fichier, sep="\t") #mettre le séparateur car ne lit pas 


        print(f"Colonnes disponibles dans {fichier} : {df.columns.tolist()}")

    except Exception as e:
        print(f" Erreur de lecture pour {fichier} : {e}")
        continue
    
    # Vérifier si la colonne INFO existe bien avant d'extraire les valeurs
    if "INFO" not in df.columns:
        print(f" Fichier {fichier} ignoré : colonne INFO manquante !")
        continue
#si une colonne essentiel info (qui contient ["pos", "len", "svtype", "svlen", "alt", "ref", "coverage", "Af"]) manque on affiche un message d'erreur et on continue
    
    condition = "hot" if "hot" in fichier else "cold" #détermination de la condition et du passage
    passage = "P50" if "P50" in fichier else "P90"
    
    if passage not in mutations:
        mutations[passage] = {}
    if condition not in mutations[passage]: #c'est ici on fait le stockage dans mutation 
        mutations[passage][condition] = []
    
    for _, row in df.iterrows(): # on parcourt chaque ligne du fichier CSV et on stocke chaque mutation sous forme de dictionnaire.

        try : 

            info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in row["INFO"].split(";") if "=" in kv} #séparation des éléments de la colonne info 
        
            variant = {
                "pos": int(row["POS"]), 
                "svtype": info_dict.get("SVTYPE", "UNKNOWN"),
                "svlen": int(info_dict.get("SVLEN", 0)),
                "alt": row["ALT"],
                "ref": row["REF"],
                "coverage": int(info_dict.get("COVERAGE", "0").split(",")[0]),  # Correction ici
                "Af": float(info_dict.get("AF", "0.0").split(",")[0])  # Correction ici
            }
            mutations[passage][condition].append(variant)

        except Exception as e:
            print(f"Erreur lors de l'extraction des données pour une ligne dans {fichier} : {e}")
            continue    

print("Lecture des fichiers terminée !")

# 3️⃣ Définition des fonctions pour comparaison

# Fonction pour calculer l'identité de séquence
def seq_identity(s1, s2):
    identity = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            identity += 1
    return identity / len(s1)

# Fonction pour comparer deux variants
def variant_equal(v1, v2, sim_thresold=1):
    if v1["svtype"] != v2["svtype"]: #On compare deux mutations pour voir si elles sont du même type. les délétions avec les délétions, les insertions avec les insertions
        return False
    
    full_length = max(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - min(v1["pos"], v2["pos"])
    common_length = max(0, min(v1["pos"] + abs(v1["svlen"]), v2["pos"] + abs(v2["svlen"])) - max(v1["pos"], v2["pos"]))
    shared = common_length / full_length if full_length > 0 else 0

    if v1["svtype"] == "INS" and shared > 0 and v1["alt"] != "<INS>" and v2["alt"] != "<INS>":
        first, second = (v1, v2) if v1["pos"] < v2["pos"] else (v2, v1)
        common_start = second["pos"] - first["pos"]
        seq1, seq2 = first["alt"][common_start:], second["alt"]
        common_stop = min(len(seq1), len(seq2))
        shared *= seq_identity(seq1[:common_stop], seq2[:common_stop])
    
    return shared >= sim_thresold

print(" Comparaison des mutations entre P50 et P90...")

# 4️⃣ Comparaison des mutations entre P50 et P90
mutations_apparues = [] #liste pour stocker les mutations suivant le status 
mutations_disparues = []
mutations_conservees = []

for condition in conditions:
    if "P50" not in mutations or "P90" not in mutations:
        print(" Erreur : Les mutations pour P50 ou P90 sont manquantes !") #Si P50 ou P90 n'existe pas dans le dictionnaire mutations, 
                                                            #cela signifie qu’on ne peut pas faire la comparaison → on arrête le script (exit()).
        exit()
    
    mutations_P50 = mutations["P50"].get(condition, [])
    mutations_P90 = mutations["P90"].get(condition, [])
    
    for variant_P90 in mutations_P90:
        found = any(variant_equal(variant_P90, v) for v in mutations_P50) #on cherche si le variant existe déjà dans P50 grâce au code elio(variant_equal)
        if found:
            mutations_conservees.append({"condition": condition, **variant_P90}) #si oui on l'ajoute à la liste mutation_conservees
        else:
            mutations_apparues.append({"condition": condition, **variant_P90}) #si non mutation apparrues à p90
    
    for variant_P50 in mutations_P50:
        found = any(variant_equal(variant_P50, v) for v in mutations_P90) #on cherche toujours grâce au code elio si la mutation à p50 existe encore à P90
        if not found:
            mutations_disparues.append({"condition": condition, **variant_P50}) #si non on l'ajoute à la liste mutations_disparues 

print(" Comparaison terminée !")

# 6️⃣ Lecture du fichier ORF.fasta et stockage des ORF
orf_data = []
fichier_orf = input("Entrez le chemin du fichier ORF.fasta : ")

# Vérifier si le chemin donné est un dossier au lieu d'un fichier
if os.path.isdir(fichier_orf):  # Vérifie si c'est un dossier
    fichiers_fasta = [f for f in os.listdir(fichier_orf) if f.endswith(".fasta")]
    if fichiers_fasta:
        fichier_orf = os.path.join(fichier_orf, fichiers_fasta[0])  # Prend le premier fichier trouvé
        print(f"Fichier ORF détecté : {fichier_orf}")
    else:
        print("Erreur : Aucun fichier .fasta trouvé dans le dossier.")
        exit()

# Vérifier si le fichier existe
if not os.path.exists(fichier_orf):
    print("Erreur : Le fichier ORF spécifié n'existe pas.")
    exit()


if os.path.exists(fichier_orf):
    for record in SeqIO.parse(fichier_orf, "fasta"):
        desc = record.description
        match = re.search(r'location=(complement\()?(\d+)\.\.(\d+)', desc)  # Récupère les positions ORF
        if match:
            orf_start = int(match.group(2))
            orf_end = int(match.group(3))
            orf_name = re.search(r'\[protein=(.+?)\]', desc)  # Récupère le nom de l'ORF
            orf_name = orf_name.group(1) if orf_name else "Unknown_ORF"
            orf_data.append({"name": orf_name, "start": orf_start, "end": orf_end})
        else:
            print(f"Erreur : Position non trouvée pour {desc}")
    print(" Lecture des ORF terminée !")
else:
    print(" Fichier ORF.fasta introuvable. Aucune vérification des ORF ne sera effectuée.")


# 7️⃣ Vérification des mutations dans les ORF

def check_orf(mutations):
    for mutation in mutations:
        mutation["ORF"] = "None"  # Valeur par défaut si la mutation n'est pas dans un ORF
        for orf in orf_data:
            if orf["start"] <= mutation["pos"] <= orf["end"]:
                mutation["ORF"] = orf["name"]  # Associer l'ORF trouvée
                break  # Arrêter dès qu'un ORF est trouvé
    return mutations


mutations_apparues = check_orf(mutations_apparues)
mutations_disparues = check_orf(mutations_disparues)
mutations_conservees = check_orf(mutations_conservees)

#  8️⃣  Sauvegarde des résultats. 

#enrégistrement en fonction de la condition

comparaison_dir = os.path.join(dossier_filtrage, "comparaison_ORF_mutation")
if not os.path.exists(comparaison_dir):
    os.makedirs(comparaison_dir)
print(f"Les résultats seront enregistrés dans {comparaison_dir}")    



print(" Mise à jour des fichiers CSV avec ORF...")
df_apparues_cold = pd.DataFrame([m for m in mutations_apparues if m["condition"] == "cold"])
df_apparues_hot = pd.DataFrame([m for m in mutations_apparues if m["condition"] == "hot"])
df_disparues_cold = pd.DataFrame([m for m in mutations_disparues if m["condition"] == "cold"])
df_disparues_hot = pd.DataFrame([m for m in mutations_disparues if m["condition"] == "hot"])
df_conservees_cold = pd.DataFrame([m for m in mutations_conservees if m["condition"] == "cold"])
df_conservees_hot = pd.DataFrame([m for m in mutations_conservees if m["condition"] == "hot"])

df_apparues_cold.to_csv(os.path.join(comparaison_dir, "mutations_apparues_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_apparues_hot.to_csv(os.path.join(comparaison_dir, "mutations_apparues_hot_ORF.csv"), index=False, sep=",", quoting=1)
df_disparues_cold.to_csv(os.path.join(comparaison_dir, "mutations_disparues_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_disparues_hot.to_csv(os.path.join(comparaison_dir, "mutations_disparues_hot_ORF.csv"), index=False, sep=",", quoting=1)
df_conservees_cold.to_csv(os.path.join(comparaison_dir, "mutations_conservees_cold_ORF.csv"), index=False, sep=",", quoting=1)
df_conservees_hot.to_csv(os.path.join(comparaison_dir, "mutations_conservees_hot_ORF.csv"), index=False, sep=",", quoting=1)



print(" ✅ Résultats mis à jour avec les ORF et enregistrés ! ")
