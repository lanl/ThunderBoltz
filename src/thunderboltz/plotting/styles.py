"""Labels, formats, plot styles, etc."""
import itertools

ls_dashed1 = (0, (1.05, 1.6, 6.5, 1.6))
ls_dashed2 = (0, (4, 2))
ls_dots1 = (0, (2.05, 1.5, 1.5, 3.05, 2.05, 1.42))
ls_dots2 = (0, (1.65, 3.5, 3.05, 2.05))
ls4 = (0, (3.5, 3.5))
ls10 = (0, (.4, 1.1, .4, 1.1, .4, 2))

# Default color cycle
DCC = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
        "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#2fb", "#00F", "#F0F",
        "#039", "#00a", "#801"]


# General default cycles
DCYC = {
    "ls": [ls_dashed1, ls_dashed2, ls_dots1, ls_dots2, ls4, ls10],
    "c": DCC,
    "marker": [".", "o", "^", "3", "1"],
}

# Prefabricated style cycles
presets = {
    "data_class": {
        "experimental": itertools.cycle({"markersize": 6, "marker": "o", "ls": "",
                        "alpha": 0.9, "c": c} for c in ["red", "green",
                        "orange", "blue", "#afa", "#aa0852", "#dba", "#abd"])},
    "name": {
        "Bolsig": itertools.cycle({"ls": "--", "c": c, "marker": None, "fillstyle": "none"} for c in DCC),
        "ThunderBoltz": itertools.cycle({"ls": "", "marker": "v", "c": c, "markersize": 12} for c in DCC),
        "Aleph": itertools.cycle([{"ls": "", "marker": "v", "c": "black", "markersize": 9}])},
}


def label_maps(s):
    """Mappings for neat math representations of parameters and column names."""
    # Non math text
    text = {
        "name": "",
        "eadf": "EADF",
        "eesd": "EESDF",
        "egen": "e$^{-}$ gen",
        "pc_scale": "Scaled Param",
        "kucukarpaci et al 1981.": 'K' + u'\u00FC' + 'c' + u'\u00FC'  + r"karpaci $et$ $al.$ 1981",
        # "el_rot": "elastic Euler rot",
        # "ion_rot": "ionizaiton Euler rot",
    }
    if s in text:
        return text[s]

    if "." in s:
        # Assume this is an experiment, remove the period at the end.
        s = s[:-1]


    if "et al" in s:
        return s.replace("et al", "$et$ $al.$")

    # Base mappings
    math_text = { # Symbol, unit symbol
        "a_n": r"\alpha/n_{\rm gas}\rm{\ (m}^{2}\rm{)}",
        "a_n_-21": r"\alpha/n_{\rm gas}\rm{\ (10}^{-21}{\rm m}^{2}\rm{)}",
        "a_n_bulk": r"\alpha/n_{\rm gas}\rm{\ \ bulk}",
        "a_n_bulk_fit": r"\alpha/N\rm{\ \ bulk}",
        "MEe": r"\langle E \rangle",
        "DT": "\Delta t",
        "DE": r"\Delta \epsilon",
        "mobN": r"\mu n_{\rm gas}\rm{\ (10}^{24}\rm{V}^{-1}\rm{m}^{-1}\rm{s}^{-1}\rm{)}",
        "mobN_bulk": r"\mu n_{\rm gas}\rm{\ (10}^{24}\rm{V}^{-1}\rm{m}^{-1}\rm{s}^{-1}\rm{)\ bulk}",
        "mobN_bulk_fit": r"\mu n_{\rm gas}\rm{\ (10}^{24}\rm{V}^{-1}\rm{m}^{-1}\rm{s}^{-1}\rm{)\ bulk}",
        "L": r"\rm{L (m)}",
        "Ered": r"\rm{E/n_{\rm gas}\ (Td)}",
        "k_tot": r"\rm{k}_{tot}",
        "name": r"\rm{}",
        "k": r"\rm{k\ (m}^3\rm{/s)}",
    }
    # Ignore if no label
    if s.replace("_std", "").replace("_err", "") not in math_text.keys():
        return s


    def alter(s):
        news, typ = "_".join(s.split("_")[:-1]), s.split("_")[-1]
        if typ == "std":
            return r"\sigma_{" + alter(news) + "}"
        elif typ == "err":
            alt = r"\%\ \rm{err}\ " + alter(news)
            # Remove units
            alt = re.sub("\(.+\)", "", alt)
            return alt
        elif s in math_text:
            # Base case, no alterations, just use label mapping
            return math_text[s]

    # Alterations such as standard deviation and error
    # return alter(s)
    return f"${alter(s)}$"

def value_maps(key):
    vmaps = {
        "DT": lambda x: f"{x:.2e}",
        "growth": lambda x: {2.0: "SST", 1.0: "PT", 4.0: "GE"}[x],
        "pc_scale": lambda x: "None" if not len(x) else str(x)
    }

    if key in vmaps:
        return vmaps[key]
    return lambda x: x
