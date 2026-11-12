# FHDeX documentation build configuration.
#
# Modelled on amrex/Docs/sphinx_documentation/source/conf.py so that the two
# documentation sets build the same way and can link to each other.

from datetime import datetime

project = "FHDeX"
author = "The FHDeX developers"
copyright = "%d, %s" % (datetime.today().year, author)
version = release = "1.0"
today_fmt = "%d %B %Y"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
]

# Cross-references into the AMReX documentation, e.g.
#   :ref:`amrex:sec:basics:multifab`
# In the PDF these become hyperlinks to the AMReX web pages.
intersphinx_mapping = {
    "amrex": ("https://amrex-codes.github.io/amrex/docs_html/", None),
    "amrex_tutorials": ("https://amrex-codes.github.io/amrex/tutorials_html/", None),
}

root_doc = "index"
source_suffix = {".rst": "restructuredtext"}
exclude_patterns = []
pygments_style = "sphinx"
numfig = True
math_number_all = True

# -- HTML ----------------------------------------------------------------------
# Built by .github/workflows/docs.yml ("make html") and deployed to the
# gh-pages branch under docs_html/.
html_theme = "sphinx_rtd_theme"
html_title = "FHDeX Documentation"
html_short_title = "FHDeX"
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
}
html_context = {
    # "Edit on GitHub" link in the top-right corner of every page
    "display_github": True,
    "github_user": "AMReX-FHD",
    "github_repo": "FHDeX",
    "github_version": "main",
    "conf_py_path": "/Docs/sphinx_documentation/source/",
}
html_static_path = ["_static"]
html_css_files = ["theme_overrides.css"]
html_show_sourcelink = False

# -- LaTeX / PDF -------------------------------------------------------------
latex_documents = [
    (root_doc, "fhdex.tex", "FHDeX User Guide", author, "manual"),
]
latex_elements = {
    "papersize": "a4paper",
    "pointsize": "10pt",
    # no index: it needs makeindex/latexmk and adds nothing for a guide this size
    "printindex": "",
}
latex_domain_indices = False
latex_show_urls = "no"
