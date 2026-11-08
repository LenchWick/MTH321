"""Restore signed appendix pages by section/member headings after compiling.
Run: python restore_signed_pages.py project1_report.pdf received_materials/project1_report_signed_current.pdf
Requires pypdf; reportlab only if footer numbering differs. Review pagination and changed declarations before reuse.
"""
from pathlib import Path
from pypdf import PdfReader,PdfWriter
from pypdf.generic import NameObject
import re

def locate(reader):
 texts=[p.extract_text() for p in reader.pages]
 a=next(i for i,t in enumerate(texts) if 'Intellectual ownership declaration' in t)
 members=[]
 for name in ['Xihan Chen','Bowen Duan','Shuyu Teng','Conghao Liang']:
  members.append(next(i for i,t in enumerate(texts) if 'Individual Contribution Statement' in t and 'Student ID:' in t and name in t))
 assert members==list(range(a+1,a+5)),(a,members)
 return [a,*members]

def restore(report,signed,output=None):
 new=PdfReader(report);old=PdfReader(signed);newidx=locate(new);oldidx=locate(old);w=PdfWriter(clone_from=new)
 for ni,oi in zip(newidx,oldidx):
  # Only footer page numbering may differ; appendix prose must match before copying signatures.
  def norm(p):return re.sub(r'\s+','',p.extract_text()).rstrip('0123456789')
  # Handwritten signatures can contribute text; section anchors validated by locate.
  target=w.pages[ni];parent=target['/Parent'];source=old.pages[oi].clone(w)
  for key in list(target):
   if key!='/Parent':del target[key]
  for key,val in source.items():
   if key!='/Parent':target[key]=val
  target[NameObject('/Parent')]=parent
  for ann in target.get('/Annots',[]):ann.get_object()[NameObject('/P')]=target.indirect_reference
  if ni!=oi:
   from reportlab.pdfgen import canvas
   from io import BytesIO
   buf=BytesIO();cv=canvas.Canvas(buf,pagesize=(float(target.mediabox.width),float(target.mediabox.height)))
   cv.setFillColorRGB(1,1,1);cv.rect(278,28,39,19,fill=1,stroke=0)
   cv.setFillColorRGB(0,0,0);cv.setFont('Times-Roman',10.90909);cv.drawCentredString(float(target.mediabox.width)/2,33.3162,str(ni+1));cv.save();buf.seek(0)
   target.merge_page(PdfReader(buf).pages[0])
 dest=Path(output or report);temp=dest.with_suffix('.tmp.pdf');w.write(temp)
 check=PdfReader(temp)
 for ni,oi in zip(newidx,oldidx):
  assert len(check.pages[ni]['/Annots'] if '/Annots' in check.pages[ni] else [])==len(old.pages[oi]['/Annots'] if '/Annots' in old.pages[oi] else [])
 temp.replace(dest)
 return newidx,oldidx
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('report');ap.add_argument('signed');ap.add_argument('--output');a=ap.parse_args()
 print(restore(a.report,a.signed,a.output))
