"""Sphinx configuration, staged by scripts/build_docs.py."""

from importlib.metadata import version as package_version

project = "GAUGE"
author = "GAUGE developers"
release = package_version("gauge_sc")
version = release
extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_design",
]
html_theme = "furo"
html_title = "GAUGE"
html_favicon = "docs/_static/brand/favicon.png"
html_static_path = ["docs/_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "light_logo": "brand/gauge_sc-logo-light.svg",
    "dark_logo": "brand/gauge_sc-logo-dark.svg",
    "sidebar_hide_name": True,
    "light_css_variables": {
        "color-brand-primary": "#245b88",
        "color-brand-content": "#245b88",
        "color-foreground-primary": "#111111",
    },
    "dark_css_variables": {
        "color-brand-primary": "#8ac0ee",
        "color-brand-content": "#8ac0ee",
    },
}
html_show_copyright = False
html_copy_source = False
exclude_patterns = [
    "docs/conf.py",
    "docs/_static/**",
    "docs/assets/README.md",
    "docs/assets/agent-audit-summary.md",
    "downloads/**",
]
myst_enable_extensions = ["colon_fence", "dollarmath", "html_admonition", "html_image"]
myst_heading_anchors = 4
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_in_temp = True
nb_execution_raise_on_error = True
nb_code_prompt_show = "Show {type}"
nb_code_prompt_hide = "Hide {type}"
autodoc_typehints = "none"
autodoc_member_order = "bysource"
napoleon_use_param = True
napoleon_use_rtype = False
