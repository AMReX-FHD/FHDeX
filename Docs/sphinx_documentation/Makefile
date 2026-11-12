# Makefile for the FHDeX Sphinx documentation (adapted from AMReX).
#
#   make html        web site in build/html       (what the GitHub Action runs)
#   make latexpdf    PDF in build/latex/fhdex.pdf  (needs a LaTeX installation)
#
# Every target first regenerates the runtime-parameter tables in
# source/params/ from the FHDeX sources (gen_params.py).

PYTHON       ?= python3
SPHINXOPTS    =
SPHINXBUILD   = $(PYTHON) -msphinx
SPHINXPROJ    = FHDeX
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help clean params Makefile

params:
	$(PYTHON) -I gen_params.py ../.. $(SOURCEDIR)/params

clean:
	rm -rf $(BUILDDIR) $(SOURCEDIR)/params

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile params
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
