#!/usr/bin/env python3
"""Generate the runtime-parameter reference for the FHDeX user guide.

FHDeX declares its input parameters in ``src_*/<name>_namespace.H`` as
``extern`` variables, sets defaults in ``Initialize<Name>Namespace()`` inside
the matching ``*_functions.cpp`` and then reads them with
``pp.query/queryarr/get/getarr("key", var)``.  This script parses those three
pieces and writes one reStructuredText table per namespace.

Usage:  gen_params.py <FHDeX root> <output dir>   (run by "make" in this directory)
"""

import re
import sys
from pathlib import Path

# (namespace label, header, source with the initializer, initializer name)
NAMESPACES = [
    ("common", "src_common/common_namespace.H",
     "src_common/common_functions.cpp", "InitializeCommonNamespace"),
    ("gmres", "src_gmres/gmres_namespace.H",
     "src_gmres/gmres_functions.cpp", "InitializeGmresNamespace"),
    ("multispec", "src_multispec/multispec_namespace.H",
     "src_multispec/multispec_functions.cpp", "InitializeMultispecNamespace"),
    ("compressible", "src_compressible/compressible_namespace.H",
     "src_compressible/compressible_functions.cpp", "InitializeCompressibleNamespace"),
    ("chemistry", "src_chemistry/chemistry_namespace.H",
     "src_chemistry/chemistry_functions.cpp", "InitializeChemistryNamespace"),
    ("reactDiff", "src_reactDiff/reactDiff_namespace.H",
     "src_reactDiff/reactDiff_functions.cpp", "InitializeReactDiffNamespace"),
    ("MFsurfchem", "src_MFsurfchem/MFsurfchem_namespace.H",
     "src_MFsurfchem/MFsurfchem_functions.cpp", "InitializeMFSurfchemNamespace"),
    ("dsmc", "src_particles/dsmc_namespace.H",
     "src_particles/dsmc_functions.cpp", "InitializeDsmcNamespace"),
    ("surfchem_mui", "src_surfchem_mui/surfchem_mui_namespace.H",
     "src_surfchem_mui/surfchem_mui_functions.cpp", "InitializeSurfChemMUINamespace"),
]

EXTERN_RE = re.compile(r"^\s*extern\s+(?P<type>.+?)\s+(?P<name>\w+)\s*;\s*(?://\s*(?P<comment>.*))?$")
QUERY_RE = re.compile(r'pp\.(?P<fn>query|queryarr|get|getarr)\(\s*"(?P<key>[^"]+)"\s*,\s*(?P<var>[\w:]+)')
ASSIGN_RE = re.compile(r"^\s*(?P<var>\w+)\s*(?P<idx>\[[^\]]*\])?\s*=\s*(?P<val>[^;=][^;]*);\s*(?://\s*(?P<comment>.*))?$")


# array sizes in the declarations -> what they mean to a user
SIZE_NAMES = {
    "AMREX_SPACEDIM": "dim",
    "3": "3",
    "MAX_SPECIES": "species",
    "MAX_SPECIES*MAX_SPECIES": "species²",
    "MAX_ELEMENT": "species pairs",
    "MAX_REACTION": "reactions",
}


def clean_type(t):
    """Turn a C++ declaration type into short user-facing notation.

    GpuArray<Real, MAX_SPECIES> -> "Real × species", Vector<int> -> "int list".
    """
    t = re.sub(r"\s+", " ", t.replace("AMREX_GPU_MANAGED", "").replace("amrex::", "")).strip()
    m = re.fullmatch(r"GpuArray<\s*(\w+)\s*,\s*([\w*]+)\s*>", t)
    if m:
        return "%s × %s" % (m.group(1), SIZE_NAMES.get(m.group(2), m.group(2)))
    m = re.fullmatch(r"Array2D<\s*(\w+)\s*,.*MAX_SPECIES.*MAX_SPECIES.*>", t)
    if m:
        return "%s × species²" % m.group(1)
    m = re.fullmatch(r"Array2D<\s*(\w+)\s*,.*MAX_REACTION.*MAX_SPECIES.*>", t)
    if m:
        return "%s × reactions × species" % m.group(1)
    m = re.fullmatch(r"Vector<\s*(\w+)\s*>", t)
    if m:
        return "%s list" % m.group(1)
    return {"IntVect": "int × dim", "std::string": "string"}.get(t, t)


def parse_header(path):
    """Return {var: (type, trailing comment, comment block above)}."""
    lines = path.read_text().splitlines()
    out = {}
    for i, line in enumerate(lines):
        m = EXTERN_RE.match(line)
        if m:
            out[m.group("name")] = (clean_type(m.group("type")),
                                    (m.group("comment") or "").strip(),
                                    block_comment(lines, i))
    return out


def function_body(src, func):
    """Return the text of ``void func(...) { ... }`` (brace matched)."""
    m = re.search(r"void\s+%s\s*\([^)]*\)\s*\{" % func, src)
    if not m:
        raise SystemExit("initializer %s not found" % func)
    depth, i = 1, m.end()
    while depth and i < len(src):
        depth += {"{": 1, "}": -1}.get(src[i], 0)
        i += 1
    return src[m.end():i - 1]


DEF_RE = re.compile(r"^\s*(?:AMREX_GPU_MANAGED\s+)?[\w:<>, ]+?\s+\w+::(?P<var>\w+)\s*(?:=\s*(?P<val>[^;]+))?;\s*(?://\s*(?P<comment>.*))?$")


def block_comment(lines, i):
    """Return the // comment block directly above lines[i].

    A block followed by more code is usually a section heading
    ("// Time-step control"), not a description of the first statement under
    it, so the block is only used when the declaration stands on its own,
    i.e. the next line is blank or a closing brace.
    """
    if i + 1 < len(lines):
        nxt = lines[i + 1].strip()
        if nxt and not nxt.startswith("}"):
            return ""
    out = []
    j = i - 1
    while j >= 0 and lines[j].strip().startswith("//"):
        text = lines[j].strip().lstrip("/").strip()
        if text:
            out.insert(0, text)
        j -= 1
    return " ".join(out)


def parse_initializer(path, func, decls):
    """Return [(key, name, fn, default, (trailing, block))] in read order.

    The default is the last plain assignment to the variable that precedes
    its pp.query in the initializer, or else the initializer on its
    file-scope definition.  Array parameters are often read into a temporary
    (``pp.queryarr("prob_lo", temp)``), so the key is used as the variable
    name whenever it matches a declaration.
    """
    src = path.read_text()
    file_lines = src.splitlines()
    file_defaults, file_comments = {}, {}
    for i, line in enumerate(file_lines):
        m = DEF_RE.match(line)
        if m:
            if m.group("val"):
                file_defaults[m.group("var")] = m.group("val").strip()
            file_comments[m.group("var")] = ((m.group("comment") or "").strip(),
                                             block_comment(file_lines, i))

    body = function_body(src, func)
    lines = body.splitlines()
    offsets, pos = [], 0
    for line in lines:
        offsets.append(pos)
        pos += len(line) + 1
    assigns = []  # (offset, var, value, comment)
    for i, line in enumerate(lines):
        m = ASSIGN_RE.match(line)
        if m and "pp." not in line:
            val = m.group("val").strip()
            if m.group("idx"):
                val += " (all)"
            c = ((m.group("comment") or "").strip(), block_comment(lines, i))
            assigns.append((offsets[i], m.group("var"), val, c))

    rows, seen = [], set()
    for m in QUERY_RE.finditer(body):
        key, fn = m.group("key"), m.group("fn")
        if key in seen:
            continue
        seen.add(key)
        name = key if key in decls else m.group("var").split("::")[-1]
        default, comments = file_defaults.get(name), file_comments.get(name, ("", ""))
        for off, var, val, c in assigns:
            if off > m.start():
                break
            if var == name:
                default, comments = val, c
        rows.append((key, name, fn, default, comments))
    return rows


def useful(comment):
    """Drop placeholder comments such as "default value" or "x (no default)"."""
    c = comment.strip().lower()
    return bool(c) and c not in ("default value", "default") and "(no default)" not in c


def rst_escape(s):
    return s.replace("*", r"\*").replace("|", r"\|").replace("`", r"\`")


def write_table(label, header, source, func, root, outdir):
    decls = parse_header(root / header)
    rows = parse_initializer(root / source, func, decls)

    lines = [
        ".. This file is generated by Docs/sphinx_documentation/gen_params.py -- do not edit.",
        "",
        ".. _params:%s:" % label,
        "",
        "``%s``" % label,
        "=" * (len(label) + 4),
        "",
        "Declared in ``%s``; defaults set in ``%s()`` (``%s``)."
        % (header, func, source),
        "",
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 36 14 12 38",
        "   :class: longtable",
        "",
        "   * - Parameter",
        "     - Type",
        "     - Default",
        "     - Description",
    ]
    for key, name, fn, default, (init_trail, init_block) in rows:
        typ, hdr_trail, hdr_block = decls.get(name, ("", "", ""))
        # prefer a comment on the same line, header first, then the comment
        # block above the declaration or default value
        comment = next((c for c in (hdr_trail, init_trail, hdr_block, init_block)
                        if useful(c)), "")
        if fn in ("get", "getarr"):
            default = "*required*"
        else:
            default = "``%s``" % default if default else "--"
        lines += [
            "   * - ``%s``" % key,
            "     - %s" % typ if typ else "     -",
            "     - %s" % default,
            "     - %s" % rst_escape(comment) if comment else "     -",
        ]
    (outdir / ("%s.rst" % label)).write_text("\n".join(lines) + "\n")
    return len(rows)


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    root, outdir = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    for ns in NAMESPACES:
        n = write_table(*ns, root=root, outdir=outdir)
        print("  %-14s %4d parameters" % (ns[0], n))


if __name__ == "__main__":
    main()
