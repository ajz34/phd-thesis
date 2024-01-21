import bibtexparser

with open("thesis.bib", "r") as f:
    lines = f.read().split("\n")

for id in range(len(lines)):
    if lines[id].strip().startswith("keywords"):
        lines[id] = ""

with open("thesis.bib", "w") as f:
    f.write("\n".join(lines))

