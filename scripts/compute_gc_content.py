"""
Calcule le GC% de chaque espece a partir d'un alignement nucleotidique
au format NEXUS (style seaview), le croise avec un fichier de
temperatures (NEXUS continuous), et ecrit un fichier NEXUS combine
a 2 colonnes (Temp, GC) pret a etre utilise dans RevBayes.

Usage :
    python3 compute_gc_content.py

Entrees attendues (a adapter si besoin) :
    - archaea.nex   : alignement ADN, noms de taxons eventuellement
                       tronques (format seaview "[n] NOM\\nSEQUENCE")
    - archaea.temp  : temperatures, format NEXUS CONTINUOUS
                       (une ligne "espece<TAB>valeur" par taxon)
    - archaea.tree  : arbre Newick, utilise seulement pour recuperer
                       les noms complets des taxons (matching avec
                       les noms tronques de archaea.nex)

Sortie :
    - archaea_temp_gc.nex : NEXUS CONTINUOUS, colonnes "Temp" et "GC"
"""

import re
import sys

# --------------------------------------------------------------
# Chemins des fichiers : passes en arguments de ligne de commande,
# sinon valeurs par defaut (cherches dans le dossier courant).
#
#   python3 compute_gc_content.py archaea.nex archaea.temp archaea.tree [sortie.nex]
# --------------------------------------------------------------
nex_path  = sys.argv[1] if len(sys.argv) > 1 else "archaea.nex"
temp_path = sys.argv[2] if len(sys.argv) > 2 else "archaea.temp"
tree_path = sys.argv[3] if len(sys.argv) > 3 else "archaea.tree"
out_path  = sys.argv[4] if len(sys.argv) > 4 else "archaea_temp___gc.nex"

# Caracteres valides dans une sequence nucleotidique (IUPAC + gaps)
_SEQ_CHARS = set("ACGTURYKMSWBDHVNacgturykmswbdhvn-?.")


def _looks_like_sequence(token):
    return len(token) > 0 and all(ch in _SEQ_CHARS for ch in token)


# --------------------------------------------------------------
# 1) Lecture des sequences ADN depuis le fichier NEXUS
#
#    Supporte automatiquement TROIS formats de bloc MATRIX :
#    (a) format seaview avec index : "[n] NOM" sur une ligne, puis
#        la sequence repliee sur plusieurs lignes suivantes ;
#    (b) format "une ligne par taxon" : "NOM<espaces/tab>SEQUENCE" ;
#    (c) format "nom seul" : le nom de l'espece seul sur sa ligne
#        (sans crochets), suivi de la sequence repliee sur plusieurs
#        lignes.
#
#    La detection se fait sans mots-cles fixes : toute ligne qui
#    n'est PAS composee uniquement de caracteres de sequence valides
#    (A/C/G/T/U/IUPAC/gap) est consideree comme un nom de taxon.
#    C'est robuste car les noms d'especes contiennent presque
#    toujours un "_" ou des minuscules hors de l'alphabet IUPAC.
# --------------------------------------------------------------
def read_nexus_dna(path):
    sequences = {}
    current = None
    with open(path) as f:
        in_matrix = False
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.strip().upper().startswith("MATRIX"):
                in_matrix = True
                continue
            if not in_matrix:
                continue
            s = line.strip()
            if s == "" or s == ";":
                continue
            if s.upper() == "END;":
                break

            # format (a) : en-tete "[n] NOM"
            m = re.match(r'^\[\d+\]\s*(\S+)\s*$', s)
            if m:
                current = m.group(1)
                sequences[current] = ""
                continue

            # format (b) : "NOM<espace(s)>SEQUENCE" sur une seule ligne
            parts = re.split(r'\s+', s, maxsplit=1)
            if len(parts) == 2 and len(parts[1]) > 30 and _looks_like_sequence(parts[1]):
                name, seq = parts
                sequences[name] = seq
                current = None  # entree complete, pas de continuation attendue
                continue

            # format (c) : ligne qui n'est pas une sequence valide
            # -> c'est un nom de taxon seul sur sa ligne
            if not _looks_like_sequence(s):
                current = s
                sequences[current] = ""
                continue

            # sinon : ligne de continuation d'une sequence
            # (formats a et c)
            if current is not None:
                sequences[current] += s

    return sequences


# --------------------------------------------------------------
# 2) Lecture des temperatures (fichier NEXUS CONTINUOUS)
# --------------------------------------------------------------
def read_temp_nexus(path):
    temps = {}
    with open(path) as f:
        in_matrix = False
        for line in f:
            line = line.strip()
            if line.upper().startswith("MATRIX"):
                in_matrix = True
                continue
            if not in_matrix or line in (";", "") or line.upper() == "END;":
                continue
            parts = line.split()
            if len(parts) == 2:
                try:
                    temps[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return temps


# --------------------------------------------------------------
# 3) Lecture des noms de taxons complets depuis l'arbre Newick
#    (utilises pour faire correspondre les noms tronques de
#    archaea.nex aux noms complets de archaea.temp)
# --------------------------------------------------------------
def read_tree_tip_names(path):
    with open(path) as f:
        tree_str = f.read()
    names = re.findall(r'[\(,]([A-Za-z][A-Za-z0-9_\.]*)\s*:', tree_str)
    return sorted(set(names))


# --------------------------------------------------------------
# 4) Correspondance nom tronque (archaea.nex) <-> nom complet
#    (archaea.temp / archaea.tree), par comparaison de prefixe
#    insensible a la casse
# --------------------------------------------------------------
def match_truncated_names(sequences, tip_names):
    def norm(s):
        return s.upper().replace("-", "_")

    mapping = {}
    for nexname in sequences:
        nn = norm(nexname).rstrip("_")
        candidates = []
        for full in tip_names:
            fn = norm(full)
            L = min(len(nn), len(fn))
            if nn[:L] == fn[:L]:
                candidates.append(full)
        if len(candidates) != 1:
            raise ValueError(f"Correspondance ambigue ou introuvable pour "
                              f"'{nexname}' : candidats = {candidates}")
        mapping[nexname] = candidates[0]
    return mapping


# --------------------------------------------------------------
# 5) Calcul du GC% d'une sequence
#    (gaps "-" et caracteres ambigus exclus du denominateur)
# --------------------------------------------------------------
def gc_content(seq):
    seq = seq.upper()
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    total = a + c + g + t
    other = len(seq) - total  # gaps, N, ?, etc.
    gc = (g + c) / total if total > 0 else None
    return gc, total, other


# --------------------------------------------------------------
# Programme principal
# --------------------------------------------------------------
def main():
    sequences = read_nexus_dna(nex_path)
    if len(sequences) == 0:
        raise RuntimeError(
            f"Aucune sequence trouvee dans '{nex_path}'. "
            f"Verifie le chemin du fichier et son format "
            f"(bloc MATRIX attendu, format seaview '[n] NOM' ou "
            f"'NOM<TAB>SEQUENCE' sur une ligne)."
        )

    temps = read_temp_nexus(temp_path)
    if len(temps) == 0:
        raise RuntimeError(f"Aucune temperature trouvee dans '{temp_path}'.")

    tip_names = read_tree_tip_names(tree_path)
    if len(tip_names) == 0:
        raise RuntimeError(f"Aucun nom de taxon trouve dans '{tree_path}'.")

    mapping = match_truncated_names(sequences, tip_names)

    rows = []
    for nexname, seq in sequences.items():
        full = mapping[nexname]
        gc, total, other = gc_content(seq)
        t = temps.get(full)
        rows.append((full, t, gc, total, other))

    rows.sort(key=lambda r: r[0])

    # ---- rapport dans la console ----
    print(f"{'Espece':40s} {'Temp':>7s} {'GC%':>8s} "
          f"{'sites_utilises':>14s} {'gaps/N':>7s}")
    for full, t, gc, total, other in rows:
        t_str = f"{t:7.2f}" if t is not None else f"{'?':>7s}"
        print(f"{full:40s} {t_str} {gc:8.3f} {total:14d} {other:7d}")

    missing_temp = [r[0] for r in rows if r[1] is None]
    if missing_temp:
        print("\nATTENTION - especes sans temperature connue :", missing_temp)

    # ---- ecriture du fichier NEXUS combine ----
    with open(out_path, "w") as out:
        out.write("#NEXUS\n\n")
        out.write("BEGIN DATA;\n")
        out.write(f"DIMENSIONS NTAX={len(rows)} NCHAR=2;\n")
        out.write("FORMAT DATATYPE=CONTINUOUS MISSING=?;\n")
        out.write("CHARSTATELABELS 1 Temp, 2 GC;\n")
        out.write("MATRIX\n")
        for full, t, gc, total, other in rows:
            t_str = f"{t:.4f}" if t is not None else "?"
            gc_str = f"{gc:.4f}" if gc is not None else "?"
            out.write(f"{full}\t{t_str}\t{gc_str}\n")
        out.write(";\nEND;\n")

    print(f"\nFichier ecrit : {out_path}")


if __name__ == "__main__":
    main()