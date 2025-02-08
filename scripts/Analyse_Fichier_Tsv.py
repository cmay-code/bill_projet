import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns #extension de Matplotlib qui permet de créer des graphiques plus lisibles et stylisés avec moins de code.
from fpdf import FPDF #générer un pdf avec toutes les stats et graphes 

# Détection automatique du chemin du script et du dossier des fichiers filtrés
script_dir = os.path.dirname(os.path.abspath(__file__))  
filter_dir = os.path.join(script_dir, "..", "data", "fusion_results", "filter_variant")

# Vérifier que le dossier existe
if not os.path.exists(filter_dir):
    print(f"Erreur : Le dossier '{filter_dir}' n'existe pas.")
    exit()


# ------------------------------------------ #
# CONFIGURATION DU DOSSIER DE SORTIE         #
# ------------------------------------------ #

output_dir = os.path.join(script_dir, "..", "data", "analyses_mutations")  # Dossier où seront enregistrés les fichiers
os.makedirs(output_dir, exist_ok=True)  # Crée le dossier s'il n'existe pas




# ------------------------------------------ #
# SCRIPT DE COMPARAISON DES MUTATIONS       #
# ------------------------------------------ #
# Script pour comparer les mutations entre les passages P50 et P90
# pour identifier celles qui ont disparu ou apparu.
# Il vérifie aussi la présence de ces mutations dans P15.
# Ce script traite à la fois les conditions FROID et CHAUD.
# ------------------------------------------ #



# Liste des conditions à traiter (FROID et CHAUD)
conditions = ["froid", "chaud"]

for condition in conditions:
    print(f" Analyse des mutations pour la condition : {condition.upper()}")
    
    # Définition des fichiers d'entrée pour chaque condition
    p50_file = os.path.join(filter_dir, f"P50_{condition}_filtered.tsv")
    p90_file = os.path.join(filter_dir, f"P90_{condition}_filtered.tsv")
    p15_file = os.path.join(filter_dir, f"P15_{condition}_filtered.tsv")

    
    # Chargement des données en tant que DataFrame
    p50 = pd.read_csv(p50_file, sep="\t")  # Lecture des mutations à P50
    p90 = pd.read_csv(p90_file, sep="\t")  # Lecture des mutations à P90
    p15 = pd.read_csv(p15_file, sep="\t")  # Lecture des mutations à P15
    
    # Sélection des colonnes importantes pour la comparaison
    cols = ["CHROM", "POS", "REF", "ALT"]
    p50 = p50[cols]  # Garder uniquement les colonnes nécessaires
    p90 = p90[cols]
    p15 = p15[cols]


    
    
    # ------------------------------------------ #
    # IDENTIFICATION DES MUTATIONS DISPARUES    #
    # ------------------------------------------ #
    # Mutations présentes à P50 mais absentes à P90
    mutations_disparues = p50.merge(p90, on=cols, how='left', indicator=True)
    mutations_disparues = mutations_disparues[mutations_disparues['_merge'] == 'left_only'].drop(columns=['_merge'])
    
    # Vérifier si ces mutations étaient présentes à P15
    mutations_disparues["Présente_P15"] = mutations_disparues.apply(
        lambda row: ((p15[cols] == row[cols]).all(axis=1)).any(), axis=1)
    
    # Sauvegarde du fichier des mutations disparues
    mutations_disparues.to_csv(os.path.join(output_dir, f"mutations_disparues_{condition}.tsv"), sep="\t", index=False)
    print(f"Fichier 'mutations_disparues_{condition}.tsv' généré dans 'data/' !")


    
    # ------------------------------------------ #
    # IDENTIFICATION DES MUTATIONS APPARUES     #
    # ------------------------------------------ #
    # Mutations présentes à P90 mais absentes à P50
    mutations_apparees = p90.merge(p50, on=cols, how='left', indicator=True)
    mutations_apparees = mutations_apparees[mutations_apparees['_merge'] == 'left_only'].drop(columns=['_merge'])
    
    # Vérifier si ces mutations étaient présentes à P15
    mutations_apparees["Présente_P15"] = mutations_apparees.apply(
        lambda row: ((p15[cols] == row[cols]).all(axis=1)).any(), axis=1)
    
    # Sauvegarde du fichier des mutations apparues
    mutations_apparees.to_csv(os.path.join(output_dir, f"mutations_apparees_{condition}.tsv"), sep="\t", index=False)
    print(f"Fichier 'mutations_apparees_{condition}.tsv' généré dans 'data/' !")
    
print(" Comparaison des mutations terminée pour toutes les conditions !")




# ------------------------------------------ #
# SCRIPT D'ANNOTATION DES MUTATIONS         #
# ------------------------------------------ #
# Script pour comparer les mutations entre les passages P50 et P90
# pour identifier celles qui ont disparu ou apparu.
# Ensuite, il annote ces mutations en fonction des ORFs
# en utilisant un fichier GFF contenant les annotations génomiques.
# ------------------------------------------ #


# Liste des types de mutations à analyser
types_mutations = ["disparues", "apparues"]

# Chargement du fichier GFF contenant les annotations génomiques
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, "..", "data")

gff_file = os.path.join(data_dir, "genomic.gff")
 # Chemin du fichier GFF

def load_gff_annotations(gff_file):
    """
    Fonction pour lire le fichier GFF et extraire les régions codantes (CDS).
    Retourne un DataFrame contenant les positions des CDS et leurs annotations.
    """
    annotations = []  # Liste pour stocker les annotations
    
    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#"):  # Ignorer les commentaires
                continue
            
            parts = line.strip().split("\t")  # Séparer les colonnes
            if len(parts) < 9:
                continue  # Vérifier que la ligne contient bien des informations complètes
            
            feature_type = parts[2]  # Type de la région (CDS, gene, etc.)
            if feature_type != "CDS":  # On ne garde que les séquences codantes
                continue
            
            chrom = parts[0]  # Chromosome
            start = int(parts[3])  # Position de début
            end = int(parts[4])  # Position de fin
            attributes = parts[8]  # Métadonnées associées
            
            gene_name = "Unknown"
            for attr in attributes.split(";"):
                if attr.startswith("gene="):
                    gene_name = attr.split("=")[1]  # Extraire le nom du gène
            
            annotations.append([chrom, start, end, gene_name])
    
    # Conversion en DataFrame
    gff_df = pd.DataFrame(annotations, columns=["CHROM", "START", "END", "GENE"])
    return gff_df

# Charger les annotations génomiques
gff_data = load_gff_annotations(gff_file)
print(" Annotations GFF chargées avec succès !")

for condition in conditions:
    for mutation_type in types_mutations:
        print(f" Annotation des mutations {mutation_type} pour la condition : {condition.upper()}")
        
        # Définition des fichiers de mutations
        mutations_file = os.path.join(output_dir, f"mutations_{mutation_type}_{condition}.tsv")
        
        # Charger les mutations
        try:
            mutations = pd.read_csv(mutations_file, sep="\t")
        except FileNotFoundError:
            print(f" Fichier {mutations_file} introuvable, passage...")
            continue
        
        # Vérification des colonnes importantes
        cols = ["CHROM", "POS", "REF", "ALT"]
        mutations = mutations[cols]  # Sélectionner les colonnes essentielles
        
        # Déterminer si chaque mutation est dans un CDS
        def annotate_mutation(row):
            """Fonction qui vérifie si une mutation est dans une séquence codante."""
            in_cds = gff_data[(gff_data["CHROM"] == row["CHROM"]) &
                              (gff_data["START"] <= row["POS"]) &
                              (gff_data["END"] >= row["POS"])].copy()
            
            if not in_cds.empty:
                return in_cds.iloc[0]["GENE"]  # Retourne le premier gène trouvé
            return "Intergénique"  # Si la mutation est en dehors des CDS
        
        # Appliquer l'annotation aux mutations
        mutations["Annotation"] = mutations.apply(annotate_mutation, axis=1)
        
        # Sauvegarde des mutations annotées
        output_file = os.path.join(output_dir, f"mutations_annotées_{mutation_type}_{condition}.tsv")
        mutations.to_csv(output_file, sep="\t", index=False)
        print(f"Fichier '{output_file}' généré !")

print(" Annotation des mutations terminée pour toutes les conditions et types de mutations !")



# ------------------------------------------ #
# SCRIPT D'ANALYSE DES MUTATIONS           #
# ------------------------------------------ #
# Script pour comparer les mutations entre les passages P50 et P90
# pour identifier celles qui ont disparu ou apparu.
# Ensuite, pour annoter ces mutations en fonction des ORFs
# et génère des visualisations pour mieux comprendre les tendances.
# ------------------------------------------ #


# Dictionnaire pour stocker les statistiques de répartition
repartition_stats = []

# ------------------------------------------ #
# LECTURE DES FICHIERS ANNOTÉS ET CALCUL DES STATISTIQUES
# ------------------------------------------ #


for condition in conditions:
    for mutation_type in types_mutations:
        file_name = f"mutations_annotées_{mutation_type}_{condition}.tsv"
        try:
            # Chargement du fichier des mutations annotées
            df = pd.read_csv(file_name, sep="\t")
            total_mutations = len(df)  # Nombre total de mutations
            
            # Compter les mutations intergéniques (hors gènes)
            intergenic_count = len(df[df["Annotation"] == "Intergénique"])
            
            # Le reste des mutations sont dans des gènes codants (CDS)
            coding_count = total_mutations - intergenic_count
            
            # Stocker les résultats dans une liste pour analyse
            repartition_stats.append({
                "Condition": condition,
                "Type": mutation_type,
                "Intergéniques": intergenic_count,
                "Codantes (CDS)": coding_count
            })
        except FileNotFoundError:
            print(f" Fichier introuvable : {file_name}")

# Conversion en DataFrame pour la visualisation
repartition_df = pd.DataFrame(repartition_stats)
# Sauvegarde du tableau de statistiques
tableau_stats = "stats_globales.tsv"
repartition_df.to_csv(tableau_stats, sep="\t", index=False)
print(f"Fichier '{tableau_stats}' généré !")




# ------------------------------------------------------------------------ #
# PLOT : RÉPARTITION DES MUTATIONS ENTRE CDS ET INTERGÉNIQUES
# ------------------------------------------------------------------------- #

plt.figure(figsize=(10, 6))  # Définir la taille de la figure

# Transformer le DataFrame pour faciliter la visualisation avec Seaborn

melted_df = repartition_df.melt(id_vars=["Condition", "Type"], 
                                value_vars=["Intergéniques", "Codantes (CDS)"],
                                var_name="Catégorie",
                                value_name="Nombre de mutations")

# Création du graphique en barres avec Seaborn
sns.barplot(x="Condition", y="Nombre de mutations", hue="Catégorie", data=melted_df, palette=["#2ecc71", "#f39c12"])

# Ajout des titres et labels
plt.xlabel("Condition (Froid vs Chaud)")  # Label axe X
plt.ylabel("Nombre de mutations")  # Label axe Y
plt.title("Répartition des mutations entre CDS et Intergéniques")  # Titre du graphique
plt.legend(title="Type de Mutation")  # Légende pour différencier CDS vs Intergénique
plt.xticks(rotation=0)  # Rotation des étiquettes X pour meilleure lisibilité
plt.grid(axis="y", linestyle="--", alpha=0.7)  # Ajout d'une grille pour lisibilité

# Sauvegarde du graphique en image
plt.savefig(os.path.join(output_dir, "repartition_mutations.png"))

# Affichage du graphique
plt.show()

print(" Graphique de répartition des mutations généré et sauvegardé !")




# ------------------------------------------ #
# ANALYSE DES GÈNES LES PLUS MUTÉS
# ------------------------------------------ #
# Création d'un dictionnaire pour stocker le comptage des mutations par gène
gene_mutation_counts = []

for condition in conditions:
    for mutation_type in types_mutations:
        file_name = os.path.join(output_dir, f"mutations_annotées_{mutation_type}_{condition}.tsv")
        try:
            # Chargement du fichier
            df = pd.read_csv(file_name, sep="\t")
            
            # Compter le nombre de mutations par gène
            gene_counts = df["Annotation"].value_counts().reset_index()
            gene_counts.columns = ["Gene", "Nombre de Mutations"]
            
            # Ajouter l'information de condition et type de mutation
            gene_counts["Condition"] = condition
            gene_counts["Type de Mutation"] = mutation_type
            
            # Ajouter aux résultats globaux
            gene_mutation_counts.append(gene_counts)
        except FileNotFoundError:
            print(f" Fichier introuvable : {file_name}")

# Fusionner les résultats en un DataFrame final
gene_mutation_df = pd.concat(gene_mutation_counts, ignore_index=True)

# Sauvegarde du fichier récapitulatif
gene_mutation_df.to_csv("genes_mutations_summary.tsv", sep="\t", index=False)
print("Fichier 'genes_mutations_summary.tsv' généré !")



# ------------------------------------------ #
# PLOT : GÈNES LES PLUS MUTÉS
# ------------------------------------------ #
plt.figure(figsize=(12, 6))

# Sélectionner les gènes les plus mutés
top_genes = gene_mutation_df.groupby("Gene")["Nombre de Mutations"].sum().nlargest(10).reset_index()

# Tracer un graphique des 10 gènes les plus mutés
sns.barplot(x="Nombre de Mutations", y="Gene", data=top_genes, palette="viridis")

# Ajouter des titres et labels
plt.xlabel("Nombre de Mutations")
plt.ylabel("Gène")
plt.title("Top 10 des gènes les plus mutés")

# Sauvegarder l'image
plt.savefig("top_genes_mutations.png")
plt.show()

print(" Graphique des gènes les plus mutés généré et sauvegardé !")


# ------------------------------------------ #
# GÉNÉRATION DU RAPPORT PDF
# ------------------------------------------ #
pdf = FPDF()
pdf.set_auto_page_break(auto=True, margin=15)
pdf.add_page()
pdf.set_font("Arial", style='B', size=16)
pdf.cell(200, 10, "Rapport d'Analyse des Mutations", ln=True, align='C')
pdf.ln(10)

# Ajout des statistiques
pdf.set_font("Arial", size=12)
pdf.cell(0, 10, "Statistiques Globales", ln=True, align='L')
pdf.ln(5)
for index, row in repartition_df.iterrows():
    pdf.cell(0, 10, f"{row['Condition']} - {row['Type']}: {row['Total Mutations']} mutations", ln=True)
pdf.ln(5)

# Ajout du graphique des gènes les plus mutés
pdf.cell(0, 10, "Graphique: Gènes les plus mutés", ln=True)
pdf.ln(5)
pdf.image("top_genes_mutations.png", x=10, w=180)
pdf.ln(10)

# Ajout du graphique de répartition des mutations
pdf.cell(0, 10, "Graphique: Répartition des mutations entre CDS et Intergéniques", ln=True)
pdf.ln(5)
pdf.image("repartition_mutations.png", x=10, w=180)
pdf.ln(10)

# SAUVEGARDE DU PDF DANS `data/`
pdf.output(os.path.join(output_dir, "Rapport_Analyse_Mutations.pdf"))
print(" Rapport PDF généré dans 'data/' !")

