"""Generate the GAUGE adaptive-canopy SVG assets with outlined Inter lettering."""

from pathlib import Path
import shutil

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs" / "_static" / "brand"


def canopy(blue, light_blue, orange, *, compact=False):
    """A shared tapered trunk with coarse crowns and finer lower branches."""
    trunk = (
        f'<path fill="{blue}" d="M82 177 C88 155 87 129 88 101 '
        'L88 49 C88 45 96 45 96 49 L96 102 C96 128 97 154 102 177Z"/>'
    )
    branches = [
        ("M91 108 C87 92 62 85 49 66", 8.5),
        ("M95 136 C112 111 107 98 115 82 C121 73 133 66 140 57", 8.5),
        ("M98 133 C114 113 130 107 150 97", 8),
    ]
    if not compact:
        branches += [
            ("M90 160 C80 137 74 125 57 123 C40 123 29 113 22 103", 7),
            ("M57 123 C51 116 47 109 47 102", 4.5),
            ("M49 122 C40 124 32 126 27 127", 4.5),
            ("M99 161 C113 145 127 147 134 131 C137 125 141 122 146 121", 5.5),
            ("M116 147 C124 147 129 145 133 141", 4.5),
        ]
    stem_paths = "\n".join(
        f'<path d="{d}" stroke-width="{width}"/>' for d, width in branches
    )
    crowns = [
        (92, 32, 28, blue),
        (39, 58, 24, blue),
        (144, 52, 21, light_blue),
        (156, 94, 15, orange),
    ]
    if not compact:
        crowns += [
            (20, 99, 11, light_blue),
            (47, 99, 7, light_blue),
            (26, 127, 7.5, light_blue),
            (146, 121, 8, light_blue),
            (132, 141, 7.5, light_blue),
        ]
    crown_paths = "\n".join(
        f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color}"/>'
        for x, y, r, color in crowns
    )
    return (
        trunk + f'\n<g fill="none" stroke="{blue}" stroke-linecap="round" '
        f'stroke-linejoin="round">\n{stem_paths}\n</g>\n{crown_paths}'
    )


def outlined_wordmark(font, text, *, x, baseline, size, color):
    glyphs = font.getGlyphSet()
    cmap = font.getBestCmap()
    units = font["head"].unitsPerEm
    paths = []
    cursor = 0
    for character in text:
        name = cmap[ord(character)]
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(pen)
        paths.append(
            f'<path d="{pen.getCommands()}" transform="translate({cursor} 0)"/>'
        )
        cursor += font["hmtx"][name][0] - 24
    return (
        f'<g fill="{color}" transform="translate({x} {baseline}) '
        f'scale({size / units} {-size / units})">\n' + "\n".join(paths) + "\n</g>",
        cursor * size / units,
    )


def svg(contents, width, *, title="GAUGE", description=""):
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.2f} 184" '
        f'role="img" aria-label="{title}">\n<title>{title}</title>\n'
        f"<desc>{description}</desc>\n{contents}\n</svg>\n"
    )


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    font = TTFont(ROOT / "src/gauge_sc/fonts/Inter-Bold.otf")
    colors = {
        "light": ("#0072B2", "#56B4E9", "#E69F00", "#111111"),
        "dark": ("#75C7EE", "#B4DFF4", "#F4B942", "#F8FAFC"),
    }
    description = (
        "An adaptive canopy: broad crowns and fine branches share a tapered trunk. "
        "One crown is highlighted in orange."
    )
    for mode, (blue, light_blue, orange, ink) in colors.items():
        mark = canopy(blue, light_blue, orange)
        lettering, word_width = outlined_wordmark(
            font, "GAUGE", x=205, baseline=135, size=108, color=ink
        )
        (OUTPUT / f"gauge_sc-logo-{mode}.svg").write_text(
            svg(
                mark + "\n" + lettering,
                205 + word_width + 20,
                description=description,
            )
        )
        (OUTPUT / f"gauge_sc-icon-{mode}.svg").write_text(
            svg(mark, 184, title="GAUGE canopy", description=description)
        )
        (OUTPUT / f"gauge_sc-favicon-{mode}.svg").write_text(
            svg(
                canopy(blue, light_blue, orange, compact=True),
                184,
                title="GAUGE canopy",
            )
        )
    shutil.copyfile(
        ROOT / "src/gauge_sc/fonts/LICENSE.txt", OUTPUT / "Inter-LICENSE.txt"
    )
    (OUTPUT / "README.txt").write_text(
        "GAUGE adaptive-canopy mark, redrawn as native SVG from the author's approved reference.\n"
        "Lettering uses Inter Bold outlines and requires no installed fonts.\n"
        "Light: blue #0072B2, sky blue #56B4E9, orange #E69F00, ink #111111.\n"
        "Dark: blue #75C7EE, sky blue #B4DFF4, orange #F4B942, ink #F8FAFC.\n"
        "Logo and icon files have transparent backgrounds. Favicon SVGs omit the fine lower branches.\n"
        "Regenerate SVGs with python scripts/generate_logo.py. Font attribution is in Inter-LICENSE.txt.\n"
    )
    print(f"Generated SVG assets: {OUTPUT}")


if __name__ == "__main__":
    main()
