import pandas as pd
import os
import sys

# Vérification des arguments
if len(sys.argv) < 2:
    print("Erreur : Aucun dossier de données spécifié !")
    print("Utilisation : python mutations_spécifiques.py /chemin/vers/dossier")
    sys.exit(1)

# Récupération du chemin du dossier de comparaison
dossier_comparaison = sys.argv[1]
dossier_sortie = os.path.join(dossier_comparaison, "variants_final")

# Création du dossier variants_final s'il n'existe pas
os.makedirs(dossier_sortie, exist_ok=True)

# Fonction pour charger un fichier CSV
def charger_fichier(nom_fichier):
    """Charge un fichier CSV et affiche un avertissement s'il manque."""
    chemin = os.path.join(dossier_comparaison, nom_fichier)
    if os.path.exists(chemin):
        return pd.read_csv(chemin)
    else:
        print(f"⚠️ Fichier manquant : {nom_fichier}")
        return pd.DataFrame()  # Retourne un DataFrame vide si le fichier n'existe pas

# Charger les fichiers des mutations conservées
mut_conserves_hot_p30_p65 = charger_fichier("mutations_conservees_hot_P30_P65.csv")
mut_conserves_hot_p65_p90 = charger_fichier("mutations_conservees_hot_P65_P90.csv")
mut_conserves_cold_p30_p65 = charger_fichier("mutations_conservees_cold_P30_P65.csv")
mut_conserves_cold_p65_p90 = charger_fichier("mutations_conservees_cold_P65_P90.csv")

# Vérifier le contenu des fichiers pour ORF
print("🔍 Vérification des fichiers d'entrée")
for nom, df in [
    ("mut_conserves_hot_p30_p65", mut_conserves_hot_p30_p65),
    ("mut_conserves_hot_p65_p90", mut_conserves_hot_p65_p90),
    ("mut_conserves_cold_p30_p65", mut_conserves_cold_p30_p65),
    ("mut_conserves_cold_p65_p90", mut_conserves_cold_p65_p90)
]:
    if not df.empty:
        print(f"{nom} - ORF value counts:")
        print(df["ORF"].value_counts(dropna=False))

# Assurer la cohérence des types
colonnes_fusion = ["pos", "ORF"]

for df in [mut_conserves_hot_p30_p65, mut_conserves_hot_p65_p90, mut_conserves_cold_p30_p65, mut_conserves_cold_p65_p90]:
    if not df.empty:
        df["pos"] = df["pos"].astype(str)
        df["ORF"] = df["ORF"].astype(str)

# Exclure les ORF "None" avant la fusion pour HOT
mut_conserves_hot_p30_p65 = mut_conserves_hot_p30_p65[mut_conserves_hot_p30_p65["ORF"] != "None"]
mut_conserves_hot_p65_p90 = mut_conserves_hot_p65_p90[mut_conserves_hot_p65_p90["ORF"] != "None"]

# Fusionner les mutations stables
mutations_stables_hot = mut_conserves_hot_p30_p65.merge(mut_conserves_hot_p65_p90, on=colonnes_fusion, how="inner")
mutations_stables_cold = mut_conserves_cold_p30_p65.merge(mut_conserves_cold_p65_p90, on=colonnes_fusion, how="inner")

# Vérifier après fusion
print("🔍 Après fusion:")
print(f"mutations_stables_hot: {mutations_stables_hot.shape}")
print(f"mutations_stables_cold: {mutations_stables_cold.shape}")

# Sauvegarde des mutations stables
mutations_stables_hot.to_csv(os.path.join(dossier_sortie, "mutations_stables_hot.csv"), index=False)
mutations_stables_cold.to_csv(os.path.join(dossier_sortie, "mutations_stables_cold.csv"), index=False)
print("✅ Mutations stables enregistrées.")

# Charger les fichiers des mutations apparues/disparues
mut_apparues_hot_p30 = charger_fichier("mutations_apparues_hot_P30.csv")
mut_apparues_cold_p30 = charger_fichier("mutations_apparues_cold_P30.csv")
mut_disparues_hot_p15 = charger_fichier("mutations_disparues_hot_P15.csv")
mut_disparues_cold_p15 = charger_fichier("mutations_disparues_cold_P15.csv")

# Identifier les mutations spécifiques au stress thermique
mutations_specifiques_hot = mutations_stables_hot[
    mutations_stables_hot["pos"].isin(mut_apparues_hot_p30["pos"]) & 
    ~mutations_stables_hot["pos"].isin(mut_disparues_hot_p15["pos"])
]

mutations_specifiques_cold = mutations_stables_cold[
    mutations_stables_cold["pos"].isin(mut_apparues_cold_p30["pos"]) & 
    ~mutations_stables_cold["pos"].isin(mut_disparues_cold_p15["pos"])
]

# Vérification après filtrage
print(f"mutations_specifiques_hot: {mutations_specifiques_hot.shape}")
print(f"mutations_specifiques_cold: {mutations_specifiques_cold.shape}")

# Sauvegarde des mutations spécifiques au stress thermique
mutations_specifiques_hot.to_csv(os.path.join(dossier_sortie, "mutations_specifiques_stress_hot.csv"), index=False)
mutations_specifiques_cold.to_csv(os.path.join(dossier_sortie, "mutations_specifiques_stress_cold.csv"), index=False)

print("✅ Mutations spécifiques au stress thermique enregistrées.")
