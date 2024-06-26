import re
import os
from os.path import join as pjoin
import sys

def add_laur(source, laur):
    """Edit Sphinx Latex to include LAUR number on title page."""
    with open(source) as f:
        s = f.read()
    m = re.search(r"\\date", s)
    i = s.count("\n", 0, m.start())
    lines = s.split("\n")
    laur_str = f"LA-UR-{laur}"
    lines[i] = lines[i][:-1] + r"\\" + laur_str + r"}"
    with open(source, "w") as f:
        f.write("\n".join(lines))

def set_title(source, title):
    """Edit Sphinx Latex to include flexible title."""
    with open(source) as f:
        s = f.read()
    m = re.search(r"\\title", s)
    i = s.count("\n", 0, m.start())
    lines = s.split("\n")
    lines[i] = r"\title{" + title + r"}"
    with open(source, "w") as f:
        f.write("\n".join(lines))

def add_affiliation(source):
    """Edit Sphinx Latex to include affiliation."""
    with open(source) as f:
        s = f.read()
    lines = s.split("\n")
    a = "Los Alamos National Laboratory, Los Alamos, NM, 87545"

    match = re.search("usepackage", s)
    pckg = s.count("\n", 0, match.start())
    if "authblk" not in lines[pckg+1]:
        lines.insert(pckg+1, r"\usepackage{authblk}")

    match = re.search("author", s)
    author = s.count("\n", 0, match.start())
    if "affil" not in lines[author+2]:
        lines.insert(author+2, f"\\affil{{{a}}}")

    with open(source, "w") as f:
        f.write("\n".join(lines))

def split_sphinx_table(latex_source, class_name, class_section, at):
    """Parse latex script and split the table for a given subsection
    at a given row of the table."""
    with open(latex_source) as f:
        lines = f.readlines()
    in_class = False
    in_section = False
    in_split_loc = False
    split_line = "NOT FOUND"
    ttype = "tabulary"
    add = "{\linewidth}[t]" # term for online the tabulary mode

    for i, l in enumerate(lines):
        # Find class location
        if "subsection" in l and class_name in l:
            in_class = True
        if not in_class: continue

        # Find table type (tabular or longtable)
        if "longtable" in l:
            ttype = "longtable"
            add = ""

        # Find section location
        if "subsection" in l and class_section in l:
            in_section = True
        if not in_section:
            continue

        # Find split location
        if (all(s in l for s in ["hyperref", "detokenize", "api", at])
                or ("sphinxcode" in l and at in l)):
            in_split_loc = True
        if not in_split_loc: continue

        if "end{tabulary}" in l:
            print("Exiting on already split table")
            # Don't split an already split method
            return
        if "sphinxhline" in l:
            split_line = i
            break

    if split_line == "NOT FOUND":
        print("NOT FOUND")
        return

    split_str = (
       "\sphinxbottomrule\n"
       f"\\end{{{ttype}}}\n"
       "%%%%% SPLIT TABLE %%%%%\n"
       "\\newpage\n"
       f"\\begin{{{ttype}}}{add}{{\X{{1}}{{2}}\X{{1}}{{2}}}}\n"
       "\sphinxtoprule\n"
       "\sphinxtableatstartofbodyhook\n"
    )
    lines[split_line] = split_str

    with open(latex_source, "w") as f:
        f.write("".join(lines))

def fix_doc_runoff():
    """Fix Sphinx Latex issues leading to tables that run off the page."""
    path = pjoin("docs", "build", "latex")
    file =  "thunderboltz.tex"
    file_path = pjoin(path, file)
    split_sphinx_table(file_path, "thunderboltz.ThunderBoltz", "Methods", "get_counts")
    split_sphinx_table(file_path, "thunderboltz.parameters.WrapParameters", "Attributes", "duration")
    split_sphinx_table(file_path, "thunderboltz.parameters.OutputParameters", "Attributes", "k\_ion")

def remove_member(src_path, class_name, class_section, member_name):
    with open(src_path) as f:
        lines = f.readlines()
    in_class = False
    in_section = False
    at_member = False
    lnum = None

    for i, l in enumerate(lines):
        # Find class location
        if "subsection" in l and class_name in l:
            in_class = True
        if not in_class: continue

        # Find section location
        if "subsection" in l and class_section in l:
            in_section = True
        if not in_section: continue

        if all(s in l for s in ["inxcode{", "sphinxupquote", member_name]):
            at_member = True
            break

    if at_member:
        new_lines = lines[:i-2] + lines[i+5:]

        with open(src_path, "w") as f:
            f.write("".join(new_lines))

def remove_members():
    path = pjoin("docs", "build", "latex")
    file = "thunderboltz.tex"
    file_path = pjoin(path, file)
    section = "thunderboltz.parameters.OutputParameters"
    remove_member(file_path, section, "Attributes", r"mobN")
    remove_member(file_path, section, "Attributes", r"mobN\_bulk")
    remove_member(file_path, section, "Attributes", r"a\_n")
    remove_member(file_path, section, "Attributes", r"a\_n\_bulk")
    remove_member(file_path, section, "Attributes", r"n\_gas")
    remove_member(file_path, section, "Attributes", r"k\_ion")
    remove_member(file_path, section, "Attributes", r"k\_1")

def extract_table(tfile, section, subsection):
    fpath = pjoin("docs", "build", "latex", "thunderboltz.tex")
    with open(fpath, "r") as f:
        lines = f.readlines()

    at_section = False
    for i, line in enumerate(lines):
        if "subsection" in line and section in line:
            at_section = True
        if not at_section: continue

        if "subsubsection" in line and subsection in line:
            break

    start = i+3
    shift = lines[start:].index("\n") - 1
    tlines = lines[start:start+shift]
    if "subsubsection" in tlines[-1]:
        tlines = tlines[:-1]
    if tlines:
        with open(tfile, "w") as f:
            f.write("".join(tlines))

def extract_tables():
    # Input Table
    tfile = pjoin("docs", "build", "latex", "tables", "tb_param_table.tex")
    extract_table(tfile, "TBParameters", "Attributes")
    tfile = pjoin("docs", "build", "latex", "tables", "out_param_table.tex")
    extract_table(tfile, "OutputParameters", "Attributes")
    tfile = pjoin("docs", "build", "latex", "tables", "particle_out_table.tex")
    extract_table(tfile, "ParticleParameters", "Attributes")

if __name__ == "__main__":
    path = pjoin("docs", "build", "latex")
    file = "thunderboltz.tex"
    source = pjoin(path, file)

    if "runoff" in sys.argv:
        fix_doc_runoff()

    if "remove-members" in sys.argv:
        remove_members()

    if "copy-table" in sys.argv:
        extract_tables()

    if "add_afil" in sys.argv:
        add_affiliation(source)

    if "--title" in sys.argv:
        ti = sys.argv.index("--title")
        set_title(source, sys.argv[ti+1])

    if "--laur" in sys.argv:
        ti = sys.argv.index("--laur")
        add_laur(source, sys.argv[ti+1])

    og = os.getcwd()
    os.chdir(path)
    os.system(f"pdflatex {file}")
    os.chdir(og)
