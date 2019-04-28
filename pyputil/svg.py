import io
from dataclasses import dataclass
import typing as tp

import cairosvg
import lxml.etree as etree
from PIL import Image

from pyputil.io import write_bytes

NSMAP = {'xlink': 'http://www.w3.org/1999/xlink'}


@dataclass
class ViewBounds:
    x: float
    y: float
    w: float
    h: float


class Svg:
    def __init__(self, view_bounds: ViewBounds,
                 background_color: tp.Optional[str] = None):
        self.bounds = view_bounds
        self.root = etree.Element("svg", nsmap=NSMAP, attrib={
            'width': '{:.1f}'.format(view_bounds.w),
            'height': '{:.1f}'.format(view_bounds.h),
            'version': '1.1',
            'viewBox': '{:.1f} {:.1f} {:.1f} {:.1f}'.format(
                view_bounds.x, view_bounds.y,
                view_bounds.w, view_bounds.h,
            ),
            'xmlns': "http://www.w3.org/2000/svg",
        })

        self._defs = etree.SubElement(self.root, 'defs')
        self._def_names = set()

        if background_color:
            etree.SubElement(self.root, 'rect', attrib={
                'id': 'background',
                'x': '{:.1f}'.format(self.bounds.x),
                'y': '{:.1f}'.format(self.bounds.y),
                'width': '{:.1f}'.format(self.bounds.w),
                'height': '{:.1f}'.format(self.bounds.h),
                'fill': background_color,
            })

    def add_def(self, tag, attrib: tp.Dict, parent=None):
        if parent is None:
            parent = self._defs
        if 'id' in attrib:
            self._def_names.add('#' + attrib['id'])
        return etree.SubElement(parent, tag, attrib=attrib)

    def use_def(self, parent, href: str, attrib=None):
        if attrib is None:
            attrib = {}

        if href not in self._def_names:
            raise KeyError("href {} not found in definitions".format(href))

        return etree.SubElement(parent, 'use', nsmap=NSMAP, attrib={
            **attrib,
            '{http://www.w3.org/1999/xlink}href': href,
        })

    def write(self, filename):
        tree = etree.ElementTree(self.root)
        tree.write(
            filename,
            encoding='utf-8',
            pretty_print=True,
            xml_declaration=True
        )

    def write_png(self, filename):
        write_bytes(filename, self.to_png_bytes())

    def to_bytes(self) -> bytes:
        return etree.tostring(
            etree.ElementTree(self.root),
            encoding='utf8',
            method='xml')

    def to_png_bytes(self) -> bytes:
        return cairosvg.svg2png(bytestring=self.to_bytes())

    def to_pil_image(self) -> Image.Image:
        return Image.open(io.BytesIO(self.to_png_bytes()))


def svgs_to_gif_bytes(svgs: [Svg], duration: int = 100, loop: int = 0):
    assert len(svgs) > 0
    images: [Image.Image] = [s.to_pil_image() for s in svgs]

    # then we just use PIL to write all frames to a byte string
    output = io.BytesIO()
    # note: un-optimized output so maybe use something like
    #   convert out.gif -coalesce -layers OptimizePlus opt.gif
    #   to reduce size
    images[0].save(
        output, format='GIF',
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop)

    return output.getvalue()


def write_svg_gif(
        filename: str,
        svgs: [Svg],
        duration: int = 100,
        loop: int = 0
):
    bs = svgs_to_gif_bytes(svgs=svgs, duration=duration, loop=loop)
    write_bytes(filename, bs)
