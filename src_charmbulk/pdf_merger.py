#!/usr/bin/env python3
"""
pdf_merger.py — PDF/PNG page extraction, collage, and multi-PDF merge-by-page.

Combines and refactors:
  - merge_pdf.py (collage helpers from cfAnRes/tools/template/)
  - collage_pdfs_by_page.ipynb (page-grouping workflow from cfAnRes/tools/notebook/charmbulk/)

Core functions
  extract_pdf_pages_as_images  — render all pages of multiple PDFs to PIL Images
  collage_png_pages_to_single  — arrange a list of PNG images into one collage PNG
  collage_pdf_pages_to_single  — arrange PDF pages into one collage page (PDF + optional PNG)
  collage_pdfs_by_page         — high-level: group same-page-number across PDFs → one collage PNG per page

Dependencies
  pip install pypdf pymupdf pillow
"""

import math
from pathlib import Path
from typing import Optional, Union

import fitz  # PyMuPDF
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation

# ─────────────────────────────────────────────────────────────────────
#  core collage functions: PDF page → PIL Image → collage PNG/PDF
# ─────────────────────────────────────────────────────────────────────
def collage_pdf_pages_to_single(
    input_pdf: Union[str, list[str]],
    output_pdf: str,
    output_png: Optional[str] = None,
    page_width: float = 1920,
    page_height: float = 1080,
    cols: Optional[int] = None,
    margin: float = 1.0,
    render_zoom: float = 2.0,
):
    """collage all pages of a PDF (or multiple PDFs) into one page, output as PDF and optionally PNG.

    Parameters
    ----------
    input_pdf : str | list[str]
        input PDF path or list of PDF paths to merge (order matters).
    output_pdf : str
        output PDF path.
    output_png : str | None
        if not None, also output a PNG image of the collage.
    page_width, page_height : float
        collage page size in points, default 1920×1080 (16:9).
    cols : int | None
        fixed number of columns; if None, auto-calculate based on 16:9 aspect ratio.
    margin : float
        margin in points to leave around each page.
    render_zoom : float
        zoom factor for rendering PDF pages to PNG (for output_png); higher = better quality but slower.

    Returns
    -------
    (output_pdf, output_png)
    """
    if isinstance(input_pdf, list):
        temp_writer = PdfWriter()
        for pdf_path in input_pdf:
            temp_reader = PdfReader(pdf_path)
            for page in temp_reader.pages:
                temp_writer.add_page(page)
        temp_path = Path(output_pdf).with_suffix(".temp.pdf")
        with open(temp_path, "wb") as f:
            temp_writer.write(f)
        input_pdf = str(temp_path)

    reader = PdfReader(input_pdf)
    n = len(reader.pages)
    if n == 0:
        raise ValueError("no pages found in input PDF(s).")

    aspect = page_width / page_height
    if cols is None:
        cols = math.ceil(math.sqrt(n * aspect))
    rows = math.ceil(n / cols)

    cell_w = page_width / cols
    cell_h = page_height / rows

    writer = PdfWriter()
    big = writer.add_blank_page(width=page_width, height=page_height)

    for i, src in enumerate(reader.pages):
        r, c = divmod(i, cols)

        sw = float(src.mediabox.width)
        sh = float(src.mediabox.height)

        avail_w = max(1.0, cell_w - 2 * margin)
        avail_h = max(1.0, cell_h - 2 * margin)

        scale = min(avail_w / sw, avail_h / sh)
        new_w = sw * scale
        new_h = sh * scale

        x_left = c * cell_w + (cell_w - new_w) / 2
        y_bottom = page_height - (r + 1) * cell_h + (cell_h - new_h) / 2

        tfm = (
            Transformation()
            .scale(sx=scale, sy=scale)
            .translate(tx=x_left, ty=y_bottom)
        )
        big.merge_transformed_page(src, tfm)

    output_pdf = str(output_pdf)
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdf, "wb") as f:
        writer.write(f)

    if output_png:
        doc = fitz.open(output_pdf)
        page = doc[0]
        mat = fitz.Matrix(render_zoom, render_zoom)
        pix = page.get_pixmap(matrix=mat, alpha=False)
        pix.save(output_png)
        doc.close()

    return output_pdf, output_png

# ─────────────────────────────────────────────────────────────────────
#  core collage functions: PNG page → collage PNG
# ─────────────────────────────────────────────────────────────────────
def collage_png_pages_to_single(
    input_png: list[str],
    output_png: str,
    page_width: float = 1920,
    page_height: float = 1080,
    cols: Optional[int] = None,
    margin: float = 2.0,
):
    """collage a list of PNG images into one PNG, arranged in a grid.

    Parameters
    ----------
    input_png : list[str]
        list of PNG file paths to collage (order matters).
    output_png : str
        output PNG file path.
    page_width, page_height : float
        collage page size in pixels, default 1920×1080 (16:9).
    cols : int | None
        fixed number of columns; if None, auto-calculate based on 16:9 aspect ratio.
    margin : float
        margin in pixels to leave around each image.

    Returns
    -------
    output_png : str
    """
    n = len(input_png)
    if n == 0:
        raise ValueError("input PNG list is empty.")

    aspect = page_width / page_height
    if cols is None:
        cols = math.ceil(math.sqrt(n * aspect))
    rows = math.ceil(n / cols)

    cell_w = page_width / cols
    cell_h = page_height / rows

    big_img = Image.new("RGB", (int(page_width), int(page_height)), (255, 255, 255))

    for i, png_path in enumerate(input_png):
        r, c = divmod(i, cols)

        img = Image.open(png_path)
        sw, sh = img.size

        avail_w = max(1.0, cell_w - 2 * margin)
        avail_h = max(1.0, cell_h - 2 * margin)

        scale = min(avail_w / sw, avail_h / sh)
        new_w = int(sw * scale)
        new_h = int(sh * scale)
        img_resized = img.resize((new_w, new_h), Image.Resampling.LANCZOS)

        x_left = int(c * cell_w + (cell_w - new_w) / 2)
        y_top = int(r * cell_h + (cell_h - new_h) / 2)

        big_img.paste(img_resized, (x_left, y_top))

    output_png = str(output_png)
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    big_img.save(output_png)
    return output_png


# ─────────────────────────────────────────────────────────────────────
#  grouping workflow: multiple PDFs → extract pages → collage by page → output PNGs
# ─────────────────────────────────────────────────────────────────────
def extract_pdf_pages_as_images(
    pdf_paths: list[str],
    dpi: int = 200,
) -> tuple[list[list[Image.Image]], int]:
    """Extract all pages from multiple PDFs and render them as PIL Images.

    Parameters
    ----------
    pdf_paths : list[str]
        list of PDF file paths to extract.
    dpi : int
        resolution

    Returns
    -------
    all_pages : list[list[Image.Image]]
        all_pages[pdf_idx][page_idx] = PIL.Image。
    max_pages : int
        maximum number of pages across all PDFs, used for grouping by page index later.
    """
    all_pages = []
    max_pages = 0

    for pdf_path in pdf_paths:
        p = Path(pdf_path)
        if not p.exists():
            raise FileNotFoundError(f"PDF file not found: {pdf_path}")

        doc = fitz.open(str(p))
        pages = []
        for page in doc:
            mat = fitz.Matrix(dpi / 72, dpi / 72)
            pix = page.get_pixmap(matrix=mat, alpha=False)
            img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
            pages.append(img)

        all_pages.append(pages)
        max_pages = max(max_pages, len(pages))
        doc.close()

    return all_pages, max_pages


def collage_pdfs_by_page(
    pdf_paths: list[str],
    output_dir: str,
    dpi: int = 400,
    page_w: int = 1620,
    page_h: int = 540,
    cols: Optional[int] = None,
    margin: float = 1.0,
    keep_temp: bool = False,
) -> list[str]:
    """
    Group multiple PDFs by page index and collage corresponding pages into single PNGs.

    Example
    -------
    3 PDFs with 5 pages each:
      page_000.png = PDF0.p1 + PDF1.p1 + PDF2.p1
      page_001.png = PDF0.p2 + PDF1.p2 + PDF2.p2
      ...

    Parameters
    ----------
    pdf_paths : list[str]
        list of PDF file paths to extract.
    output_dir : str
        output directory.
    dpi : int
        resolution for extracting pages from PDF.
    page_w, page_h : int
        width and height of the output collage page in pixels.
    cols : int | None
        fixed number of columns in the collage; if None, auto-calculate based on 16:9 aspect ratio.
    margin : float
        margin in pixels to leave around each page in the collage.
    keep_temp : bool
        if True, keep the temporary PNG files of individual pages; if False, delete them after creating the collage.

    Returns
    -------
    output_paths : list[str]
        list of generated PNG file paths.
    """
    all_pages, max_pages = extract_pdf_pages_as_images(pdf_paths, dpi=dpi)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    temp_dir = out / "_temp_pages"

    output_paths = []

    for page_idx in range(max_pages):
        # gather the same page index from all PDFs (if a PDF has fewer pages, skip it for that page index)
        group_pngs = []
        for pdf_idx, pages in enumerate(all_pages):
            if page_idx >= len(pages):
                continue
            temp_path = temp_dir / f"pdf{pdf_idx:02d}_pg{page_idx:03d}.png"
            temp_path.parent.mkdir(parents=True, exist_ok=True)
            pages[page_idx].save(temp_path)
            group_pngs.append(str(temp_path))

        if not group_pngs:
            continue

        output_png = out / f"page_{page_idx:03d}.png"
        collage_png_pages_to_single(
            input_png=group_pngs,
            output_png=str(output_png),
            page_width=page_w,
            page_height=page_h,
            cols=cols,
            margin=margin,
        )
        output_paths.append(str(output_png))

        # cleanup temp PNGs for this page
        if not keep_temp:
            for p in group_pngs:
                Path(p).unlink()

    # cleanup temp directory if empty
    if not keep_temp and temp_dir.exists():
        try:
            temp_dir.rmdir()
        except OSError:
            pass

    return output_paths
