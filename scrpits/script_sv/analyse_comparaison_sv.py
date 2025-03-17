import os
import pandas as pd
import sys

def charger_fichier(chemin, sep):
    """
    Charge un fichier CSV avec le bon séparateur et extrait les mutations.
    """
    try:
        df = pd.read_csv(chemin, sep=sep)
        df.columns = df.columns.str.lower()
        
        # Extraction des informations utiles (SVTYPE, SVLEN) si présentes dans INFO
        if 'info' in df.columns:
            df["svtype"] = df["info"].str.extract(r'SVTYPE=([^;]+)')
            df["svlen"] = df["info"].str.extract(r'SVLEN=([^;]+)').astype(float)
        
        return df
    except Exception as e:
        print(f"❌ Erreur lors du chargement du fichier {chemin}: {e}")
        return None

def comparer_mutations(df1, df2, df3, tolerance=5):
    """
    Conserve les mutations présentes dans les trois fichiers en tolérant de légères variations.
    """
    mutations_communes = []
    
    for _, row in df1.iterrows():
        pos, svtype, svlen = row["pos"], row.get("svtype", None), row.get("svlen", None)
        
        # Vérifier si cette mutation existe dans P65 et P90 avec une tolérance sur la position
        match_df2 = df2[(df2["pos"].between(pos - tolerance, pos + tolerance)) & (df2["svtype"] == svtype)]
        match_df3 = df3[(df3["pos"].between(pos - tolerance, pos + tolerance)) & (df3["svtype"] == svtype)]
        
        # Si la mutation est trouvée dans les trois fichiers, elle est conservée
        if not match_df2.empty and not match_df3.empty:
            mutations_communes.append(row)
    
    return pd.DataFrame(mutations_communes)

def extraire_mutations_uniques(df1, df2, df3, tolerance=5):
    """
    Conserve les mutations qui sont uniquement présentes dans P30 et absentes dans P65 et P90.
    """
    mutations_uniques = []
    
    for _, row in df1.iterrows():
        pos, svtype = row["pos"], row.get("svtype", None)
        
        match_df2 = df2[(df2["pos"].between(pos - tolerance, pos + tolerance)) & (df2["svtype"] == svtype)]
        match_df3 = df3[(df3["pos"].between(pos - tolerance, pos + tolerance)) & (df3["svtype"] == svtype)]
        
        # Si la mutation N'EXISTE PAS dans P65 et P90, elle est considérée comme unique
        if match_df2.empty and match_df3.empty:
            mutations_uniques.append(row)
    
    return pd.DataFrame(mutations_uniques)

def main(dossier_annotes, dossier_filtre):
    fichiers = {
        "p30_hot": (os.path.join(dossier_annotes, "mutations_annotées_P30_hot.csv"), ","),
        "p30_cold": (os.path.join(dossier_annotes, "mutations_annotées_P30_cold.csv"), ","),
        "p65_hot": (os.path.join(dossier_filtre, "P65_hot_filtre.csv"), "\t"),
        "p65_cold": (os.path.join(dossier_filtre, "P65_cold_filtre.csv"), "\t"),
        "p90_hot": (os.path.join(dossier_filtre, "P90_hot_filtre.csv"), "\t"),
        "p90_cold": (os.path.join(dossier_filtre, "P90_cold_filtre.csv"), "\t"),
    }
    
    # Charger tous les fichiers
    mutations = {key: charger_fichier(fname, sep) for key, (fname, sep) in fichiers.items()}
    
    # Vérification des fichiers chargés
    if any(df is None for df in mutations.values()):
        print("❌ Certains fichiers n'ont pas été correctement chargés. Vérifiez les chemins et formats.")
        return
    
    # Comparer les mutations entre P30, P65 et P90 pour HOT et COLD
    mutations_conservées_hot = comparer_mutations(mutations["p30_hot"], mutations["p65_hot"], mutations["p90_hot"])
    mutations_conservées_cold = comparer_mutations(mutations["p30_cold"], mutations["p65_cold"], mutations["p90_cold"])
    
    # Extraire les mutations uniques à P30
    mutations_uniques_hot = extraire_mutations_uniques(mutations["p30_hot"], mutations["p65_hot"], mutations["p90_hot"])
    mutations_uniques_cold = extraire_mutations_uniques(mutations["p30_cold"], mutations["p65_cold"], mutations["p90_cold"])
    
    # Création des dossiers de sortie
    dossier_sortie = os.path.join("mutation_P3_ORF", "mutations_conservées")
    os.makedirs(dossier_sortie, exist_ok=True)
    dossier_uniques = os.path.join("mutation_P3_ORF", "mutations_uniques")
    os.makedirs(dossier_uniques, exist_ok=True)
    
    # Sauvegarde des résultats dans les dossiers spécifiés
    mutations_conservées_hot.to_csv(os.path.join(dossier_sortie, "mutations_conservées_hot.csv"), index=False)
    mutations_conservées_cold.to_csv(os.path.join(dossier_sortie, "mutations_conservées_cold.csv"), index=False)
    mutations_uniques_hot.to_csv(os.path.join(dossier_uniques, "mutations_uniques_P30_hot.csv"), index=False)
    mutations_uniques_cold.to_csv(os.path.join(dossier_uniques, "mutations_uniques_P30_cold.csv"), index=False)
    
    print(f"✅ Fichiers générés dans {dossier_sortie} et {dossier_uniques} !")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("❌ Erreur : Vous devez fournir deux dossiers en argument !")
        print("Utilisation : python analyse_mutations.py <dossier_annotes> <dossier_filtre>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
