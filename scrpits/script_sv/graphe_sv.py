import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os
import sys

# 📌 Vérification des arguments
if len(sys.argv) != 2:
    print("❌ Utilisation: python graphe_sv.py <chemin_du_dossier>")
    sys.exit(1)

dossier = sys.argv[1]

# 📌 Vérification que le dossier existe
if not os.path.isdir(dossier):
    print(f"❌ Erreur : Le dossier {dossier} n'existe pas.")
    sys.exit(1)

# 📌 Détection automatique des fichiers
hot_files = {
    "P30": os.path.join(dossier, "P30_hot_filtre.csv"),
    "P65": os.path.join(dossier, "P65_hot_filtre.csv"),
    "P90": os.path.join(dossier, "P90_hot_filtre.csv")
}

cold_files = {
    "P30": os.path.join(dossier, "P30_cold_filtre.csv"),
    "P65": os.path.join(dossier, "P65_cold_filtre.csv"),
    "P90": os.path.join(dossier, "P90_cold_filtre.csv")
}

hot_common_file = os.path.join(dossier, "mutations_conservées_hot.csv")
cold_common_file = os.path.join(dossier, "mutations_conservées_cold.csv")

# 📌 Fonction pour charger les mutations avec détection du séparateur
def load_mutations(file_path):
    """ Charge les mutations d'un fichier CSV avec détection automatique du séparateur et de la bonne colonne."""
    try:
        with open(file_path, "r") as f:
            first_line = f.readline()
            sep = "," if "," in first_line else "\t"
        
        df = pd.read_csv(file_path, sep=sep)
        df.columns = df.columns.str.strip().str.lower()  # Normalisation des colonnes
        
        if "pos" in df.columns:
            mutations = set(df["pos"].astype(str))
            print(f"✅ {file_path} : {len(mutations)} mutations chargées.")
            return mutations
        
        print(f"⚠️ {file_path} : Aucune colonne 'pos' trouvée, fichier ignoré.")
        return set()
    except Exception as e:
        print(f"❌ Erreur lors de la lecture du fichier {file_path} : {e}")
        return set()

def plot_venn(set1, set2, set3, common, labels, title):
    """ Génère un diagramme de Venn avec uniquement les pourcentages et correction des exclusions."""
    plt.figure(figsize=(6,6))
    venn = venn3([set1, set2, set3], set_labels=labels)
    
    total_mutations = len(set1 | set2 | set3 | common)
    
    for subset_id in ['100', '010', '001', '110', '101', '011', '111']:
        label = venn.get_label_by_id(subset_id)
        if label:
            count = len(set1 & set2 & set3) if subset_id == '111' else int(label.get_text()) if label.get_text() else 0
            percentage = (count / total_mutations * 100) if total_mutations > 0 else 0
            label.set_text(f"{percentage:.1f}%")
            label.set_fontsize(12)
    
    plt.title(title, fontsize=14)
    plt.show()

# 📌 Chargement des mutations individuelles
hot_sets = {key: load_mutations(hot_files[key]) for key in hot_files}
cold_sets = {key: load_mutations(cold_files[key]) for key in cold_files}

# 📌 Chargement des mutations communes
hot_common = load_mutations(hot_common_file)
cold_common = load_mutations(cold_common_file)

# 📌 Correction pour calculer les mutations propres
for key in hot_sets:
    hot_sets[key] -= hot_common
for key in cold_sets:
    cold_sets[key] -= cold_common

# Génération des diagrammes
if all(len(hot_sets[key]) > 0 for key in hot_sets):
    print("📈 Génération du diagramme HOT...")
    plot_venn(hot_sets["P30"], hot_sets["P65"], hot_sets["P90"], hot_common, labels=("P30", "P65", "P90"), title="Diagramme de Venn - HOT")
else:
    print("❌ Impossible de générer le diagramme HOT : Un des ensembles est vide.")

if all(len(cold_sets[key]) > 0 for key in cold_sets):
    print("📈 Génération du diagramme COLD...")
    plot_venn(cold_sets["P30"], cold_sets["P65"], cold_sets["P90"], cold_common, labels=("P30", "P65", "P90"), title="Diagramme de Venn - COLD")
else:
    print("❌ Impossible de générer le diagramme COLD : Un des ensembles est vide.")
