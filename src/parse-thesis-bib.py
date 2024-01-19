import bibtexparser

with open("thesis.bib", "r") as f:
    bib = bibtexparser.load(f)

def try_pop(d, k):
    try:
        return d.pop(k)
    except KeyError:
        return

for key, dic in bib.entries_dict.items():
    try:
        dic["author"] = dic["author"].upper().replace(" AND ", " and ")
    except:
        continue
    if dic["ENTRYTYPE"] not in ["article", "book", "incollection"]:
        continue
    try_pop(dic, "url")
    try_pop(dic, "urldate")
    try_pop(dic, "file")
    try_pop(dic, "note")
    try_pop(dic, "abstract")

with open("thesis.bib", "w") as f:
    bibtexparser.dump(bib, f)

with open("thesis.bib", "r") as f:
    lines = f.readlines()
token = "".join([l for l in lines if not l.startswith("@comment")])
with open("thesis.bib", "w") as f:
    f.write(token)

