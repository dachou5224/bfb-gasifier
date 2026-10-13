#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从指定 Hamel 原 PDF 生成 T2 本地语料和导航索引，不实现或解释方程。

依赖：pdfplumber、pypdf；渲染时另需 Poppler 的 pdftoppm。
页序从 1 开始。此扫描版有双页、重复页和重复式号，不能用统一偏移。
原扫描不含印刷页 152–153；这两页用 --supplement 单独编为 PDF 114，不插入已有页序。
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import unicodedata
from pathlib import Path

SOURCE_SHA256 = "1e1a83a1e87c7de45f20811fe3322aa4d036c541d1c45a63cb845ddcb30990ef"
SUPPLEMENT_SHA256 = "d99d192003cd95e2382eace333eb3308a18703ea533c3f5e962b483fdfa3f8ac"
SUPPLEMENT_PDF_PAGE = 114

# 对照书目原图核实的姓名 OCR 修正，不依据外部文献或旧转录。
OCR_SURNAMES = {
    "BL": "BI", "GLL": "GIL", "HILLLGARDT": "HILLIGARDT", "JLANG": "JIANG",
    "LEPPALAHTL": "LEPPÄLAHTI", "LLLLEDAHL": "LILIEDAHL", "LLPPENS": "LIPPENS",
    "MLLLOLL": "MILIOLI", "PILLAL": "PILLAI", "PLLLAI": "PILLAI",
    "SLMELL": "SIMELL", "SLMINSKI": "SIMINSKI", "SLNGH": "SINGH",
    "SLT": "SIT", "WLRSUM": "WIRSUM",
}

# 从原文目录和正文标题整理；不从旧实现或旧转录推导。
SECTION_ROWS = """
1|1|Einleitung
1.1|1|Stand der Vergasungstechnologie
1.2|5|Ziel der Arbeit und Methodik
1.3|6|Literaturübersicht und Stand des Wissens
2|10|Modellbeschreibung
2.1|11|Aufbau und Struktur des Simulationsprogramms
2.2|15|Formulierung der Bilanzgleichungen des Zellenmodells
2.2.1|16|Bilanzierung gasförmiger Stoffströme
2.2.2|17|Bilanzierung fester Stoffströme
2.2.3|19|Bilanzierung der Energieströme
2.3|20|Berechnungsverfahren zur Lösung der Bilanzgleichungen
3|21|Modellierung der Strömung in blasenbildenden Wirbelschichtfeuerungs- und -Vergasungsreaktoren
3.1|21|Strömungsmechanik in der Wirbelschicht
3.1.1|21|Minimale Fluidisierung
3.1.2|23|Strömungsmechanik der blasenbildenden Wirbelschicht
3.1.3|30|Blasenbildende Wirbelschichten unter erhöhtem Druck
3.1.4|34|Stoffaustausch zwischen Blasen- und Suspensionsphase
3.2|36|Strömungsmechanik im Freiraum
3.2.1|36|Berechnung des Feststoffaustrags in den Freiraum
3.2.2|40|Berechnung der Partikelstartgeschwindigkeit
3.2.3|42|Gasbewegung im Freiraum
3.2.4|43|Berechnung der Partikelbewegung im Freiraum
4|50|Trocknung und Pyrolyse des Brennstoffs
4.1|51|Modellierung der Trocknung und Pyrolyse
4.1.1|54|Trocknung und Pyrolyse feuchter Braunkohle
4.1.2|56|Pyrolyse von Holz
4.2|58|Berechnung der Zusammensetzung der Pyrolyseprodukte
4.2.1|60|Koks- und Teerausbeute bei der Pyrolyse fossiler Brennstoffe
4.2.1.1|60|Zusammenstellung verschiedener Kohlen
4.2.1.2|63|Rheinische Braunkohle
4.2.1.3|66|Torf
4.2.2|67|Pyrolyseprodukte verschiedener Hölzer
4.3|69|Berechnung der Zusammensetzung der gasförmigen Pyrolyseprodukte
5|71|Modellierung der chemischen Reaktionen
5.1|72|Berechnung des Koksumsatzes
5.1.1|74|Theoretische Untersuchung zum Einfluss des Partikelalters
5.1.1.1|76|Formulierung der Bilanzgleichungen mit Partikelaltersverteilung
5.1.1.2|78|Beschreibung des reduzierten Simulationsmodells
5.1.1.3|80|Ergebnisse zur Simulation der Partikelaltersverteilung
5.1.2|86|Koksverbrennung
5.1.3|88|Heterogene Wassergasreaktion
5.1.4|89|Boudouard-Reaktion
5.1.5|89|Hydrierende Vergasung
5.2|90|Berechnung der homogenen Reaktionen in der Gasphase
5.2.1|90|Kohlenmonoxid-Oxidation
5.2.2|91|Wasserstoff-Oxidation
5.2.3|91|Methan-Oxidation
5.2.4|91|Homogene Wassergasreaktion
5.2.5|94|Methan-Reformierung
5.2.6|95|Umsetzung der Teere
6|98|Experimentelle Untersuchungen zur Pyrolyse, Verbrennung und Vergasung verschiedener Einsatzstoffe in der Wirbelschicht
6.1|98|Beschreibung der Versuchsanlage
6.1.1|98|Aufbau der Laborwirbelschicht und messtechnische Ausstattung
6.1.2|99|Charakterisierung der eingesetzten Brennstoffe
6.1.3|100|Messdatenaufbereitung
6.1.4|102|Versuchsdurchführung
6.2|105|Pyrolyseuntersuchungen in der Wirbelschicht
6.3|108|Untersuchungen der Verbrennung in der Wirbelschicht
6.3.1|108|Ermittlung globaler Verbrennungsraten in der Wirbelschicht
6.3.2|111|Untersuchung der Kinetik der Restkoksverbrennung in der Wirbelschicht
6.4|118|Untersuchung der Kinetik der Kohlendioxidvergasung in der Wirbelschicht
7|124|Simulation verschiedener Wirbelschichtvergasungsreaktoren
7.1|124|Vergleich von Messungen an ausgeführten Wirbelschichtvergasern mit Simulationsergebnissen
7.1.1|124|Erläuterung der Ergebnisdarstellung
7.1.2|125|Simulation der Hochtemperatur-Winkler-Vergasung
7.1.3|132|Simulation des Druckwirbelschichtvergasers des VTT Finnland
7.1.4|138|Simulation des atmosphärischen Wirbelschichtvergasers des VTT Finnland
7.1.5|142|Simulation des Wirbelschichtvergasers WSV 400
7.2|147|Parameterstudie zum Einfluss verschiedener Betriebsvariablen auf die Wirbelschichtvergasung
8|154|Diskussion der Ergebnisse und Ausblick
8.1|154|Beurteilung des Vergleichs zwischen Messung und Simulation
8.2|160|Weitere Modellentwicklung
A1|163|Berechnung des chemischen Gleichgewichts in der Gasphase
A2|169|Zusammenstellung der Pyrolyse-, Verbrennungs- und Vergasungsversuche
A2.1|169|Zusammenstellung der zur Berechnung der globalen Verbrennungsraten herangezogenen Versuchsfälle
A2.2|171|Zusammenstellung der zur Berechnung der Pyrolyseproduktausbeuten und der kinetischen Parameter der Restkoksverbrennung herangezogenen Versuchsfälle
A2.3|173|Zusammenstellung der zur Berechnung der kinetischen Parameter der Restkoksvergasung mit CO2 herangezogenen Versuchsfälle
A2.4|175|Grafische Darstellung der Gaskonzentrationen und der Wirbelschichttemperatur ausgewählter Versuche
"""

# 逐式定位清单。3.49、3.50 在原书 S.35、36 各出现一次，原样保留。
# 格式为：章号、印刷页号、该页实际存在的编号（不是依据连续编号补缺）。
EQUATION_ROWS = """
2|15|1-2
2|16|3-5
2|17|6
2|19|7-8
2|20|9
3|21|1-2
3|22|3-5
3|23|6
3|24|7-10
3|25|11-12
3|26|13-18
3|27|19-26
3|28|27-33
3|29|34-40
3|31|41-43
3|32|44
3|34|45-48
3|35|49-50
3|36|49-51
3|37|52-55
3|38|56-60
3|39|61-64
3|40|65-69
3|41|70-72
3|42|73
3|43|74-80
3|44|81-84
3|45|85-87
3|46|88-91
3|47|92-93
3|48|94-98
3|49|99-105
4|52|1-4
4|53|5-11
4|54|12
4|60|13
4|61|14-18
4|62|19-20
4|63|21
4|64|22-24
4|67|25
4|69|26-27
4|70|28-30
5|74|1-2
5|75|3
5|76|4-5
5|83|6-7
5|86|8-14
5|87|15-22
5|88|23-25
5|89|26-35
5|90|36-37
5|91|38-41
5|92|42-44
5|94|45-50
5|95|51-57
5|96|58-61
6|101|1
6|108|2
6|110|3
6|111|4-5
6|112|6-11
6|113|12-16
6|119|17-20
6|121|21-23
A|163|1-3
A|164|4-11
A|165|12-20
A|166|21-26
A|167|27-34
A|168|35-40
"""

# 少量可直接目视核对的简短公式，纯导航；其余保留空串而非发布损坏的 OCR 公式。
SHORT_FORMS = {
    "2.3": r"N_{Ges}=N_Z\cdot N_{Bil}",
    "3.6": r"\alpha_b=u_{b,i}/u_{mf}",
    "3.39": r"h_0=h_A+d_{b,0}/2",
    "3.40": r"h_0=d_{b,0}/2",
    "3.46": r"K_{bd}=(1/K_{bc}+1/K_{cd})^{-1}",
    "3.74": r"u_{gb}(h)=u_0+(u_{gb,0}-u_0)e^{-\beta_A h}",
    "3.75": r"\beta_A=1/(3d_{b,WS})",
    "3.83": r"u_r=u_P-u_g",
    "3.85": r"c_W=24/Re+z",
    "3.92": r"\Delta=4AC-B^2",
    "4.13": r"\mathrm{Kohle}\longrightarrow\mathrm{Metaplast}\longrightarrow\mathrm{Teer}",
    "4.28": r"G=\sum_{j=1}^{N_S}\mu_j n_j=\mathrm{Min!}",
    "4.29": r"dG=\sum_{j=1}^{N_S}(\partial G/\partial n_j)dn_j=0",
    "4.30": r"\sum_{j=1}^{N_S}a_{ij}n_j-b_i=0",
}


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def extract(pdf, pages_dir):
    import pdfplumber
    from pypdf import PdfReader

    reader = PdfReader(pdf)
    result = []
    with pdfplumber.open(pdf) as document:
        for number, page in enumerate(document.pages, 1):
            x0, y0, x1, y1 = page.bbox
            mid = (x0 + x1) / 2
            sides = []
            if page.rotation == 90:
                split = (y0 + y1) / 2
                halves = [("top", (x0, y0, x1, split)), ("bottom", (x0, split, x1, y1))]
            else:
                halves = [("left", (x0, y0, mid, y1)), ("right", (mid, y0, x1, y1))]
            for side, bounds in halves:
                cropped = page.crop(bounds)
                sides.append({"side": side,
                              "text": cropped.extract_text(x_tolerance=1, y_tolerance=3) or "",
                              "lines": cropped.extract_text_lines(x_tolerance=1, y_tolerance=3, return_chars=False)})
            raw = reader.pages[number - 1].extract_text() or ""
            if not raw.strip() and number not in (1, 2):
                raise ValueError("PDF 文本层无法提取：第 %d 页；按 T2 停止，需检查 OCR" % number)
            (pages_dir / ("%03d.txt" % number)).write_text(raw + "\n", encoding="utf-8")
            result.append({"pdf_page": number, "bbox": list(page.bbox), "rotation": page.rotation, "sides": sides})
            if number % 20 == 0:
                print("已提取 %d 页" % number, flush=True)
    write_json(pages_dir / "extraction.json", result)
    return result


def render(pdf, pages_dir, count):
    if all((pages_dir / ("%03d.png" % n)).exists() for n in range(1, count + 1)):
        return
    # 只复用本目录已有的同一原文渲染结果，不读取归档材料。
    generated = [pages_dir / ("page-%03d.png" % n) for n in range(1, count + 1)]
    if not all(p.exists() for p in generated):
        renderer = shutil.which("pdftoppm")
        if not renderer:
            raise RuntimeError("缺少 pdftoppm；请安装 Poppler 或使用 bundled runtime")
        subprocess.run([renderer, "-r", "300", "-png", str(pdf), str(pages_dir / "page")],
                       check=True, timeout=7200)
    for n, path in enumerate(generated, 1):
        path.rename(pages_dir / ("%03d.png" % n))


def printed_pair(number):
    if number <= 9:
        return {3: (None, "III"), 4: (None, "V"), 5: ("VI", "VII"),
                6: ("VIII", "IX"), 7: ("X", "XI"), 8: ("XII", "XIII"),
                9: ("XIV", "XV")}.get(number, (None, None))
    if number >= 112:
        return None, None
    if number == 10:
        return None, "1"
    start = 2 * number - 20
    if number in (85, 86):
        start -= 2
    return str(start), str(start + 1)


def location(pdf_page, side, printed_page):
    return {"pdf_page": pdf_page, "printed_page": printed_page, "side": side}


def section_for(printed_page, sections, y=None):
    if printed_page is None or not printed_page.isdecimal():
        return None
    candidates = [(int(s["printed_page_start"]), s.get("heading_y", -1), k)
                  for k, s in sections.items() if int(s["printed_page_start"]) <= int(printed_page)]
    if y is not None:
        candidates = [v for v in candidates if v[0] < int(printed_page) or v[1] <= y]
    return max(candidates, default=(0, 0, None))[2]


def make_sections(pages, reverse):
    sections = {}
    for row in SECTION_ROWS.strip().splitlines():
        key, printed, title = row.split("|", 2)
        loc = reverse[printed][0]
        lines = pages[loc["pdf_page"] - 1]["sides"][0 if loc["side"] == "left" else 1]["lines"]
        # 目录确定起始页；正文匹配用于页内标题定位，匹配不到就明确只定位到半页。
        # OCR 有 2,23、4.2.1.3Torf 等情况；只接受行首完整编号，排除页眉。
        matches = []
        y0 = pages[loc["pdf_page"] - 1]["bbox"][1]
        for line in lines:
            marker = re.match(r"^([Aa]\.?\s*[1lI2](?:[.,]\s*[0-9]+)*|[0-9]+(?:[.,]\s*[0-9]+)*)", line["text"])
            marker_key = normalized(marker[1].replace("l", "1").replace("I", "1")) if marker else None
            if marker_key == normalized(key) and line["top"] > y0 + 25:
                matches.append(line)
            elif key == "7.1.3" and line["text"].startswith("1A3 Simulation"):
                matches.append(line)  # PDF 76 标题编号的 OCR 误读，按目录和原图定位。
        entry = {"title_de": title, "pdf_page_start": loc["pdf_page"],
                 "printed_page_start": printed, "side": loc["side"],
                 "source": "原文目录 PDF 4–5；正文起始页"}
        if matches:
            entry["heading_y"] = round(matches[0]["top"], 2)
        else:
            entry["heading_y"] = -1
            entry["heading_note"] = "正文标题 OCR 不完整；起始页按原文目录定位"
        sections[key] = entry
    return sections


def make_equations(pages, reverse, sections):
    equations = {}
    for row in EQUATION_ROWS.strip().splitlines():
        chapter, printed, numbers = row.split("|")
        low, high = (numbers.split("-") + [numbers])[:2] if "-" in numbers else (numbers, numbers)
        loc = reverse[printed][0]
        lines = pages[loc["pdf_page"] - 1]["sides"][0 if loc["side"] == "left" else 1]["lines"]
        for number in range(int(low), int(high) + 1):
            key = "%s.%d" % (chapter, number)
            pat = re.compile(r"[<(]\s*" + chapter + r"\s*[.\-’, ]\s*" + str(number) + r"\s*[)>]")
            matches = [line for line in lines if pat.search(line["text"])
                       and not re.search(r"Gleichung|Bild|Tabelle|\bG[lL]\.", line["text"])]
            y = matches[0]["top"] if matches else None
            entry = dict(loc, section=section_for(printed, sections, y),
                         label_de="Gleichung (%s)" % key, short_form=SHORT_FORMS.get(key, ""),
                         transcription_status="navigation_only" if key in SHORT_FORMS else "not_transcribed",
                         locator="text_label" if matches else "page_image",
                         status="indexed")
            if matches:
                line = matches[0]
                entry["label_bbox_pt"] = [round(line[k], 2) for k in ("x0", "top", "x1", "bottom")]
            if key in equations:
                equations[key].setdefault("additional_locations", []).append(entry)
                equations[key]["ambiguous_number"] = True
                equations[key]["note"] = "原文印刷式号重复；须同时查看 S.35 和 S.36，不得仅按式号取式"
            else:
                equations[key] = entry
    return equations


def caption(lines, offset):
    collected = [lines[offset]["text"]]
    if collected[0].rstrip(" _").endswith("."):
        return join_lines(collected)
    last = lines[offset]
    for line in lines[offset + 1:offset + 8]:
        if (line["top"] - last["bottom"] > 6 or re.match(r"^(Bild|Tabelle|\d+(?:\.\d+)*\s)", line["text"])
                or line["text"].startswith(("Parameter", "Wert", "Anlagenparameter", "Betriebsparameter",
                                            "Exp.-Nr.", "Bezeichnung", "Kohle", "VTT-Druckwirbelschichtvergaser"))):
            break
        collected.append(line["text"])
        last = line
        if line["text"].endswith("."):
            break
    return join_lines(collected)


def join_lines(lines):
    text = "\n".join(lines)
    text = re.sub(r"(?<=\w)-\n(?=[a-zäöü])", "", text)
    return re.sub(r"\s+", " ", text).strip()


def make_captions(pages, sections):
    groups = {"figures": {}, "tables": {}}
    pat = re.compile(r"^(Bild|Tabelle)\s+([1-8A]\s*\.\s*\d+(?:\.\d+)*)\s*(?:\(Fortsetzung\))?\s*[:;]")
    for page in pages[9:111]:
        for side, printed in zip(page["sides"], printed_pair(page["pdf_page"])):
            for i, line in enumerate(side["lines"]):
                match = pat.search(line["text"])
                if not match:
                    continue
                key = re.sub(r"\s+", "", match[2])
                ocr_key = key
                if page["pdf_page"] == 96 and side["side"] == "right" and key == "A.23.1":
                    key = "A.2.3.1"  # 目视原图：OCR 漏了 2 和 3 之间的点。
                group = groups["figures" if match[1] == "Bild" else "tables"]
                entry = dict(location(page["pdf_page"], side["side"], printed),
                             section=section_for(printed, sections, line["top"]),
                             caption_de=caption(side["lines"], i), status="indexed",
                             caption_source="pdf_text_layer",
                             label_bbox_pt=[round(line[k], 2) for k in ("x0", "top", "x1", "bottom")])
                if ocr_key != key:
                    entry["ocr_label"] = ocr_key
                    entry["caption_de"] = entry["caption_de"].replace(ocr_key, key, 1)
                    entry["caption_source"] = "pdf_text_layer_with_image_label_correction"
                if key in group:
                    group[key].setdefault("additional_locations", []).append(entry)
                    if "Fortsetzung" not in line["text"]:
                        group[key]["ambiguous_number"] = True
                        group[key]["note"] = "原书出现相同编号的不同图表；保留所有位置，未自行重编号"
                else:
                    group[key] = entry
    for key in ("7.20", "7.21"):
        groups["figures"][key] = {"pdf_page": None, "printed_page": None, "section": "7.2",
                                  "caption_de": None, "status": "source_page_missing",
                                  "note": "正文 S.148–149 提及；扫描件缺 S.152–153，不能恢复图或图注",
                                  "mentioned_at": [location(84, "left" if key == "7.20" else "right",
                                                            "148" if key == "7.20" else "149")]}
    # 表 7.2 的续表标记不是标准图注，但数据明确跨到下一印刷页。
    groups["tables"]["7.2"].setdefault("additional_locations", []).append(
        dict(location(77, "left", "134"), caption_de="Fortsetzung Tabelle 7.2", status="indexed"))
    return groups


def normalized(text):
    text = unicodedata.normalize("NFKD", text).upper().replace("ß", "SS")
    return re.sub(r"[^A-Z0-9]", "", text).replace("1", "I").replace("0", "O")


def author_positions(text, names):
    """姓名必须占完整词，不把 Ko 匹配到 Koks、So 匹配到 sowie。"""
    tokens = list(re.finditer(r"[^\W_]+(?:[-'][^\W_]+)*", text))
    positions = []
    for i, token in enumerate(tokens):
        for width in range(1, min(5, len(tokens) - i + 1)):
            last = tokens[i + width - 1]
            if normalized(text[token.start():last.end()]) in names:
                positions.append((token.start(), last.end()))
    return positions


def is_reference_start(text):
    if text.startswith("BARON, T.;"):
        return True  # 原文 OCR 此条以分号而非冒号分隔作者和题名。
    prefix = text.split(":")[0]
    letters = [c for c in prefix if c.isalpha()]
    return (":" in text and 2 < len(prefix) < 160 and bool(letters)
            and sum(c.isupper() for c in letters) / len(letters) > .6)


def make_references(pages, sections):
    records = []
    current = None
    for page in pages[103:111]:
        for side, printed in zip(page["sides"], printed_pair(page["pdf_page"])):
            for line in side["lines"]:
                text = line["text"]
                if "Literaturverzeichnis" in text or text.isdecimal() or text == "1":
                    continue
                if is_reference_start(text):
                    current = {"lines": [], "location": location(page["pdf_page"], side["side"], printed)}
                    records.append(current)
                if current is not None:
                    current["lines"].append(text)
    references = {}
    for serial, record in enumerate(records, 1):
        full = join_lines(record["lines"])
        authors = full.split(":", 1)[0] if not full.startswith("BARON,") else full.split("The Effect", 1)[0]
        # 主键保留原书作者+年份；同姓同年不同作品另加序号，不合并文献。
        surname = re.sub(r"^1\s*", "", authors.split(",", 1)[0]).strip()
        surname = re.sub(r"\s+19\d{2}$", "", surname)
        if surname.startswith("RHEINBRAUN"):
            surname = "RHEINBRAUN AG"
        years = list(dict.fromkeys(re.findall(r"\b(?:19\d{2}|200[01])[ab]?\b", full)))
        year = years[0] if years else "undated"
        canonical = OCR_SURNAMES.get(normalized(surname), surname)
        if normalized(surname) == "WERTHER" and not years:
            # S.202 原图显示 1977，文本层把 1 误识别为反斜线。
            full = full.replace("(\\977)", "(1977)")
            years = ["1977"]
            year = "1977"
        base = normalized(canonical) + "_" + year
        key = base
        suffix = 2
        while key in references:
            key = base + "_%d" % suffix
            suffix += 1
        references[key] = {"full_citation": full, "bibliography_location": record["location"],
                           "first_author_ocr": surname, "first_author": canonical,
                           "year_tokens": years, "cited_for": [],
                           "citation_locations": [], "status": "bibliography_transcribed_from_text_layer",
                           "bibliography_order": serial}
    for ref in references.values():
        surname = normalized(ref["first_author"])
        # 正文 OCR 的同一人名变体；只影响导航匹配，不改原始书目转录。
        aliases = {"WERTHER": ["WEITHER"], "BOTTERJLL": ["BOTTERILL"],
                   "DESOUZASANTOS": ["SOUZASANTOS"]}.get(surname, [])
        if normalized(ref["first_author_ocr"]) != surname:
            aliases.append(normalized(ref["first_author_ocr"]))
        for page in pages[9:103]:
            if page["pdf_page"] == 85:
                continue  # S.148–149 的重复扫描不重复计数正文引用。
            for side, printed in zip(page["sides"], printed_pair(page["pdf_page"])):
                lines = side["lines"]
                for i, line in enumerate(lines):
                    window = join_lines([l["text"] for l in lines[i:i+3]])
                    window = re.sub(r"(?<=[A-Za-z])etal\.", " et al.", window)
                    plain = normalized(window)
                    if not any(name in plain for name in [surname] + aliases):
                        continue
                    positions = author_positions(window, [surname] + aliases)
                    positions = [(start, end) for start, end in positions
                                 if start < len(line["text"]) and (window[start].isupper()
                                     or surname.startswith(("VAN", "DE")))]
                    if not positions:
                        continue
                    if not any(re.search(r"\b" + re.escape(year[:4]) + r"[ab]?\b", window[end:end+100], re.I)
                               for _, end in positions for year in ref["year_tokens"]):
                        continue
                    entry = dict(location(page["pdf_page"], side["side"], printed),
                                 section=section_for(printed, sections, line["top"]),
                                 match_method="author_year_text_window_candidate",
                                 line_y=round(line["top"], 2))
                    if entry not in ref["citation_locations"]:
                        ref["citation_locations"].append(entry)
        ref["cited_for"] = sorted({c["section"] for c in ref["citation_locations"] if c["section"]})
        ref["matching_note"] = "仅表示作者和年份在正文邻近文本中出现；不自动证明该文献是某方程的依据。OCR、同姓同年与年份差异须按页面复核。"
        if not ref["citation_locations"]:
            ref["matching_note"] = "未自动定位正文引用；不等于原文未引用。按书目页和本地全文人工检索。"
    return references


def build_index(pages, pdf):
    page_map = {}
    reverse = {}
    for page in pages:
        pair = printed_pair(page["pdf_page"])
        sides = [side["side"] for side in page["sides"]]
        page_map[str(page["pdf_page"])] = {side: printed for side, printed in zip(sides, pair)}
        for side, printed in zip(sides, pair):
            if printed:
                reverse.setdefault(printed, []).append(location(page["pdf_page"], side, printed))
    sections = make_sections(pages, reverse)
    equations = make_equations(pages, reverse, sections)
    captions = make_captions(pages, sections)
    references = make_references(pages, sections)
    index = {"_meta": {"source": "Hamel, VDI Fortschritt-Berichte Reihe 6 Nr.469, 2001",
                        "source_filename": pdf.name, "source_sha256": SOURCE_SHA256,
                        "pdf_page_count": len(pages), "pdf_page_numbering": "1-based",
                        "render_dpi": 300, "local_pages": "refs/hamel/pages/NNN.png, NNN.txt",
                        "page_map_schema": "常规扫描页映射 left/right；PDF 98–103 为 top/bottom；null 表示无印刷页码；重复页在 reverse map 保留全部位置",
                        "coordinates": "label_bbox_pt、heading_y、line_y 为相对页面图像左上角的 point 坐标（已减去 pdfplumber 页 bbox 原点）；像素 = pt × render_dpi / 72",
                        "index_scope": "方程定位清单覆盖第2–6章和附录A1，第1/7/8章无带编号方程；全书小节、可提取图表和参考文献同时建立导航",
                        "short_form_policy": "只有少量目视核对的短式；空串表示未转录，绝不使用 OCR 公式作为科学依据",
                        "acceptance_status": "总司令抽查 20/20 通过（2026-10-13）；源缺页 S.152–153（§7.2 参数研究，非验证工况）已裁定接受"},
             "page_map": page_map, "printed_page_map": reverse, "equations": equations,
             "tables": captions["tables"], "figures": captions["figures"],
             "sections": sections, "cited_references": references,
             "source_issues": [
                 {"kind": "missing_pages", "printed_pages": ["152", "153"],
                  "affected_figures": ["7.20", "7.21"], "evidence": "PDF 86 为 S.150–151，PDF 87 为 S.154–155；原文 S.148–149 提及两图"},
                 {"kind": "duplicate_scan", "pdf_pages": [84, 85], "printed_pages": ["148", "149"]},
                 {"kind": "duplicate_equation_numbers", "numbers": ["3.49", "3.50"], "printed_pages": ["35", "36"]},
                 {"kind": "duplicate_table_numbers", "locations": {"5.1": ["74", "79"], "A.2.1.2": ["169", "170"]},
                  "note": "原书同号的不同表；索引保留所有位置，不自行重编号"},
                 {"kind": "rotated_scans", "pdf_pages": list(range(98, 104)),
                  "note": "附录试验曲线页有90度旋转，PDF显示为上下双页，使用top/bottom定位；全文保持原文本层"}]}
    index["coverage"] = {}
    for chapter in ("2", "3", "4", "5", "6", "7", "A"):
        eqs = [v for k, v in equations.items() if k.startswith(chapter + ".")]
        figs = [v for k, v in captions["figures"].items() if k.startswith(chapter + ".")]
        index["coverage"][chapter] = {
            "equation_numbers": len(eqs),
            "equation_occurrences": sum(1 + len(v.get("additional_locations", [])) for v in eqs),
            "tables": sum(k.startswith(chapter + ".") for k in captions["tables"]),
            "located_figures": sum(v["status"] == "indexed" for v in figs),
            "missing_figures": sum(v["status"] != "indexed" for v in figs)}
    return index


def to_image_coords(index, pages):
    """把 pdfplumber 坐标平移到页面图像坐标系：减去各页 bbox 的左上角原点。"""
    origin = {page["pdf_page"]: (page["bbox"][0], page["bbox"][1]) for page in pages}

    def walk(node):
        if isinstance(node, list):
            for item in node:
                walk(item)
            return
        if not isinstance(node, dict):
            return
        number = node.get("pdf_page", node.get("pdf_page_start"))
        if number in origin:
            ox, oy = origin[number]
            if node.get("label_bbox_pt"):
                x0, y0, x1, y1 = node["label_bbox_pt"]
                node["label_bbox_pt"] = [round(x0 - ox, 2), round(y0 - oy, 2), round(x1 - ox, 2), round(y1 - oy, 2)]
            for key in ("heading_y", "line_y"):
                if node.get(key) is not None:
                    node[key] = round(node[key] - oy, 2)
        for value in node.values():
            walk(value)

    walk(index)
    return index


def validate(index):
    located_7 = 21 if index["figures"]["7.20"]["status"] == "indexed" else 19
    expected = {"2": (9, 2, 5), "3": (105, 1, 8), "4": (30, 4, 13), "7": (0, 5, located_7)}
    for chapter, (equations, tables, figures) in expected.items():
        actual = index["coverage"][chapter]
        if (actual["equation_numbers"], actual["tables"], actual["located_figures"]) != (equations, tables, figures):
            raise ValueError("索引数量与逐页定位清单不符：章节 %s %s" % (chapter, actual))
    if len(index["cited_references"]) != 217:
        raise ValueError("参考文献分条与原文清单不符")
    for chapter, count in [("5", 61), ("6", 23), ("A", 40)]:
        if index["coverage"][chapter]["equation_numbers"] != count:
            raise ValueError("补充章节的式号定位不完整：%s" % chapter)
    for group in ("equations", "tables", "figures"):
        for key, entry in index[group].items():
            if entry["status"] != "indexed":
                continue
            mapped = index["page_map"][str(entry["pdf_page"])][entry["side"]]
            if mapped != entry["printed_page"]:
                raise ValueError("页码错配：%s %s" % (group, key))
    if not index["equations"]["3.49"].get("ambiguous_number"):
        raise ValueError("丢失原文重复式号")


def _insert_before(mapping, key, value, before):
    rebuilt = {}
    placed = False
    for existing, item in mapping.items():
        if existing == before and not placed:
            rebuilt[key] = value
            placed = True
        rebuilt[existing] = item
    if not placed:
        rebuilt[key] = value
    return rebuilt


def _image_box(line, origin):
    ox, oy = origin
    return [round(line["x0"] - ox, 2), round(line["top"] - oy, 2),
            round(line["x1"] - ox, 2), round(line["bottom"] - oy, 2)]


def render_supplement_png(pdf, png):
    renderer = shutil.which("pdftoppm")
    if renderer:
        subprocess.run([renderer, "-r", "300", "-png", "-f", "1", "-l", "1", str(pdf),
                        str(png.parent / "supplement")], check=True, timeout=120)
        (png.parent / "supplement-1.png").replace(png)
        return "pdftoppm 300 dpi"
    import fitz
    document = fitz.open(str(pdf))
    pixmap = document[0].get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
    pixmap.save(str(png))
    document.close()
    return "PyMuPDF %s，300 dpi；本机没有 pdftoppm" % fitz.version[0]


def place_supplement_png(pdf, png, skip_render):
    if skip_render:
        if png.exists():
            return "未渲染，沿用已有 114.png"
        return "未渲染，本地没有 114.png"
    return render_supplement_png(pdf, png)


def ingest_supplement_page(pdf, pages_dir, skip_render=False):
    """把补页写入本地 114.txt，并挂到 extraction.json。不改 1–113。skip_render 时不创建或改写 114.png。"""
    import pdfplumber
    from pypdf import PdfReader

    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if digest != SUPPLEMENT_SHA256:
        raise ValueError("补页 PDF SHA256 不同")
    reader = PdfReader(str(pdf))
    if len(reader.pages) != 1:
        raise ValueError("补页 PDF 应为单页对开扫描")
    with pdfplumber.open(str(pdf)) as document:
        page = document.pages[0]
        x0, y0, x1, y1 = page.bbox
        mid = (x0 + x1) / 2
        sides = []
        for side, bounds in (("left", (x0, y0, mid, y1)), ("right", (mid, y0, x1, y1))):
            cropped = page.crop(bounds)
            sides.append({"side": side,
                          "text": cropped.extract_text(x_tolerance=1, y_tolerance=3) or "",
                          "lines": cropped.extract_text_lines(x_tolerance=1, y_tolerance=3, return_chars=False)})
    left_header = sides[0]["lines"][0]["text"]
    right_header = sides[1]["lines"][0]["text"]
    if not left_header.startswith("152 ") or "153" not in right_header:
        raise ValueError("补页页眉不是 S.152 / S.153：%s | %s" % (left_header, right_header))
    label_20 = next(line for line in sides[0]["lines"] if line["text"].startswith("Bild 7.20"))
    label_21 = next(line for line in sides[1]["lines"] if line["text"].startswith("Bild 7.21"))
    if _image_box(label_20, (x0, y0))[:2] != [9.73, 495.76]:
        raise ValueError("Bild 7.20 标签位置与已核对的补页不符")
    if _image_box(label_21, (x0, y0)) != [421.22, 500.84, 743.23, 516.57]:
        raise ValueError("Bild 7.21 标签位置与已核对的补页不符")
    band = [line for line in sides[0]["lines"] if 478 < line["top"] < 490]
    box_20 = [round(min(line["x0"] for line in band) - x0, 2),
              round(min(line["top"] for line in band) - y0, 2),
              round(max(line["x1"] for line in band) - x0, 2),
              round(max(line["bottom"] for line in band) - y0, 2)]
    if box_20 != [9.73, 495.76, 335.08, 512.46]:
        raise ValueError("Bild 7.20 首行范围与已核对的补页不符")
    pages_dir.mkdir(parents=True, exist_ok=True)
    raw = reader.pages[0].extract_text() or ""
    (pages_dir / "114.txt").write_text(raw + "\n", encoding="utf-8")
    render_note = place_supplement_png(pdf, pages_dir / "114.png", skip_render)
    record = {"pdf_page": SUPPLEMENT_PDF_PAGE, "bbox": [x0, y0, x1, y1], "rotation": 0, "sides": sides}
    cache = pages_dir / "extraction.json"
    pages = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else []
    pages = [item for item in pages if item["pdf_page"] != SUPPLEMENT_PDF_PAGE]
    pages.append(record)
    write_json(cache, pages)
    return {"filename": pdf.name, "render_note": render_note, "box_20": box_20,
            "box_21": _image_box(label_21, (x0, y0))}


def merge_supplement(index, pdf, pages_dir, skip_render=False):
    """把已核对的 S.152–153 写入索引。图注以页面图像为准，文本层漏了 ü。"""
    placed = ingest_supplement_page(pdf, pages_dir, skip_render=skip_render)
    captions = {
        "7.20": ("152", "left", placed["box_20"],
                 "Bild 7.20: Simulationsergebnisse für den Modellreaktor unter Variation des "
                 "Feuchtegehalts des zugeführten Brennstoffs.",
                 location(84, "left", "148")),
        "7.21": ("153", "right", placed["box_21"],
                 "Bild 7.21: Simulationsergebnisse für den Modellreaktor unter Variation der "
                 "Aufteilung der zugeführten Luft in Primär- und Sekundärluft.",
                 location(84, "right", "149")),
    }
    for key, (printed, side, box, caption_de, mentioned) in captions.items():
        index["figures"][key] = dict(location(SUPPLEMENT_PDF_PAGE, side, printed),
                                     section="7.2", caption_de=caption_de, status="indexed",
                                     caption_source="page_image", label_bbox_pt=box,
                                     note="图注对照页面图像；文本层漏写 ü，且 7.20 被拆成单字符行",
                                     mentioned_at=[mentioned])
    index["page_map"] = _insert_before(
        index["page_map"], "114", {"left": "152", "right": "153"}, "113")
    # _insert_before 把新键放在 113 之前；114 应紧跟 113。
    page_map = index["page_map"]
    entry_114 = page_map.pop("114")
    index["page_map"] = {}
    for key, value in page_map.items():
        index["page_map"][key] = value
        if key == "113":
            index["page_map"]["114"] = entry_114
    reverse = index["printed_page_map"]
    reverse.pop("152", None)
    reverse.pop("153", None)
    reverse = _insert_before(reverse, "152", [location(SUPPLEMENT_PDF_PAGE, "left", "152")], "154")
    index["printed_page_map"] = _insert_before(
        reverse, "153", [location(SUPPLEMENT_PDF_PAGE, "right", "153")], "154")
    index["source_issues"] = [
        {"kind": "supplement_scan", "pdf_pages": [SUPPLEMENT_PDF_PAGE],
         "printed_pages": ["152", "153"], "affected_figures": ["7.20", "7.21"],
         "sha256": SUPPLEMENT_SHA256,
         "evidence": "原扫描 PDF 86 为 S.150–151，PDF 87 为 S.154–155。补页单独编为 PDF 114，左 S.152、右 S.153。"},
        *[item for item in index["source_issues"]
          if item["kind"] not in ("missing_pages", "supplement_scan")]
    ]
    index["coverage"]["7"]["located_figures"] = 21
    index["coverage"]["7"]["missing_figures"] = 0
    meta = index["_meta"]
    meta["acceptance_status"] = ("总司令抽查 20/20 通过（2026-10-13）。"
                                 "S.152–153 于 2026-10-14 由补页纳入 PDF 114，待总司令复核。")
    meta["page_map_schema"] = (meta["page_map_schema"].split("。PDF 114")[0]
                               + "。PDF 114 为补页，左 152、右 153")
    meta["supplement"] = {"filename": placed["filename"], "sha256": SUPPLEMENT_SHA256,
                          "pdf_page": SUPPLEMENT_PDF_PAGE, "printed_pages": ["152", "153"],
                          "render": placed["render_note"],
                          "note": "原扫描 113 页不含这两页。补页不插入 PDF 86 与 87 之间，已有页序不变。PDF 与页面图像不入库。"}
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", nargs="?", type=Path, help="原 PDF 路径，不复制或提交 PDF")
    parser.add_argument("--output", type=Path, default=Path("refs/hamel"))
    parser.add_argument("--reuse-extraction", action="store_true", help="复用本次 PDF 提取的本地缓存")
    parser.add_argument("--skip-render", action="store_true", help="只重建索引；不代表图像交付已完成")
    parser.add_argument("--supplement", type=Path, help="S.152–153 补页 PDF。不替换 113 页原扫描，也不入库")
    parser.add_argument("--supplement-only", action="store_true",
                        help="不重建 1–113，只把补页并入现有索引和本地页")
    args = parser.parse_args()
    if args.supplement_only:
        if args.supplement is None:
            raise ValueError("补页需要 --supplement")
        index_path = args.output / "index.json"
        index = merge_supplement(json.loads(index_path.read_text(encoding="utf-8")),
                                 args.supplement, args.output / "pages",
                                 skip_render=args.skip_render)
        validate(index)
        write_json(index_path, index)
        print("S.152–153 已并入 PDF 114。页面图像只在本地，PDF 不入库。")
        return
    if args.pdf is None:
        parser.error("需要原 PDF 路径")
    existing = args.output / "index.json"
    if existing.exists() and args.supplement is None:
        old = json.loads(existing.read_text(encoding="utf-8"))
        if "152" in old.get("printed_page_map", {}):
            raise ValueError("现有索引已含 S.152–153。重跑须带 --supplement，否则会把补页写回缺口。")
    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError("PDF SHA256 不同：此扫描版的人工定位清单不能套用到另一版本")
    pages_dir = args.output / "pages"
    pages_dir.mkdir(parents=True, exist_ok=True)
    cache = pages_dir / "extraction.json"
    if args.reuse_extraction:
        pages = [item for item in json.loads(cache.read_text(encoding="utf-8"))
                 if item["pdf_page"] <= 113]
    else:
        pages = extract(args.pdf, pages_dir)
    if len(pages) != 113:
        raise ValueError("PDF 页数不符")
    if not args.skip_render:
        render(args.pdf, pages_dir, len(pages))
    index = to_image_coords(build_index(pages, args.pdf), pages)
    if args.supplement is not None:
        index = merge_supplement(index, args.supplement, pages_dir, skip_render=args.skip_render)
    validate(index)
    write_json(args.output / "index.json", index)
    print(json.dumps(index["coverage"], ensure_ascii=False, indent=2))
    if args.supplement is None:
        print("原 PDF 缺 S.152–153（§7.2，非验证工况）：缺口已记录并经总司令裁定接受。")
    else:
        print("S.152–153 已并入 PDF 114。页面图像只在本地，PDF 不入库。")


if __name__ == "__main__":
    main()
