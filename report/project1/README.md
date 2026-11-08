# Project 1 report and presentation

- [Latest signed report](project1_report.pdf): 39 pages, including the team AI log and four individual contribution statements. Corresponding prose revisions were integrated on 12 October 2026.
- [Slides](project1_slides.pdf): approved 17-page presentation.
- [Report LaTeX](project1_report.tex) and [slides LaTeX](project1_slides.tex).
- [Speaker script](../../notes/project1/speaker_script.md).

The PDFs are the existing final deliverables, copied byte for byte. LaTeX report resource paths are adjusted for this repository layout; its body, equations and results are unchanged. Compile from this directory with Tectonic or a compatible LaTeX distribution. Slides use `presentation_assets/` via their graphic search path.

A fresh LaTeX compile produces an unsigned report. To restore the existing signed appendix pages, keep a separate copy of the signed report, compile, and run:

```sh
python -m pip install pypdf reportlab
python restore_signed_pages.py project1_report.pdf /path/to/saved_signed_report.pdf
```

The helper locates the A4 declaration and members' ICS pages by headings. Review the appendix before reusing signatures after edits. The checked-in PDF already contains the preserved signatures; no restoration is necessary to read it.
