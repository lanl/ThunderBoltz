# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ThunderBoltz"
copyright = "2023, Ryan Park"
author = 'Ryan Park, Brett Scheiner, Mark Zammit'
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # View the documentation compilation times.
    "sphinx.ext.duration",
    # Actually run documented interpreter script
    "sphinx.ext.doctest",
    # Auto generate docs for methods, module, classes etc.
    "sphinx.ext.autodoc",
    # Auto generate summaries for classes and their corresponding
    # docs.
    "sphinx.ext.autosummary",
    # Mathtext ":math:" in docs
    "sphinx.ext.mathjax",
    # Google docstring parsing
    "sphinx.ext.napoleon",
    # Link source code in docs
    "sphinx.ext.viewcode",
    # Link to external docs
    "sphinx.ext.intersphinx",
    # RTD html theme
    "sphinx_rtd_theme",
]

intersphinx_mapping = {
    "python": ("http://docs.python.org/", None),
    "pandas": ("http://pandas.pydata.org/pandas-docs/dev", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

templates_path = ["_templates"]
exclude_patterns = []

# autoclass_content = "class"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
# html_static_path = ["_static"]
# html_css_files = ["custom.css"]

# -- Latex options
latex_elements = {
  'extraclassoptions': 'openany,oneside'
}
