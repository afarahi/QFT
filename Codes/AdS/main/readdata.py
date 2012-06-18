from numpy             import sqrt
#import easy to use xml parser called minidom:
from xml.dom.minidom import parseString
#all these imports are standard on most modern python implementations

def read_data(file_name = 'input.xml'):
   #open the xml file for reading:
   file = open(file_name,'r')
   #convert to string:
   data = file.read()
   #close file because we dont need it anymore:
   file.close()
   #parse the xml you got from the file
   dom = parseString(data)
   #retrieve the first xml tag (<tag>data</tag>) that the parser finds with name tagName:
   xmlTag = dom.getElementsByTagName('rs')[0].toxml()
   #strip off the tag (<tag>data</tag>  --->   data):
   rs  = float(xmlTag.replace('<rs>','').replace('</rs>',''))

   xmlTag = dom.getElementsByTagName('r0')[0].toxml()
   r0  = float(xmlTag.replace('<r0>','').replace('</r0>',''))

   xmlTag = dom.getElementsByTagName('R')[0].toxml()
   R   = float(xmlTag.replace('<R>','').replace('</R>',''))
   R2  = R**2

   xmlTag = dom.getElementsByTagName('kx')[0].toxml()
   kx  = float(xmlTag.replace('<kx>','').replace('</kx>',''))

   xmlTag = dom.getElementsByTagName('kt_real')[0].toxml()
   xmlTagp= dom.getElementsByTagName('kt_complex')[0].toxml()
   kt  = float(xmlTag.replace('<kt_real>','').replace('</kt_real>','')) + float(xmlTagp.replace('<kt_complex>','').replace('</kt_complex>',''))*1.0j

   xmlTag = dom.getElementsByTagName('muq_real')[0].toxml()
   xmlTagp= dom.getElementsByTagName('muq_complex')[0].toxml()
   muq = float(xmlTag.replace('<muq_real>','').replace('</muq_real>','')) + float(xmlTagp.replace('<muq_complex>','').replace('</muq_complex>',''))*1.0j

   xmlTag = dom.getElementsByTagName('Rm2')[0].toxml()
   Rm2 = float(xmlTag.replace('<Rm2>','').replace('</Rm2>',''))

   xmlTag = dom.getElementsByTagName('dimention')[0].toxml()
   d   = float(xmlTag.replace('<dimention>','').replace('</dimention>',''))

   xmlTag = dom.getElementsByTagName('h')[0].toxml()
   h   = float(xmlTag.replace('<h>','').replace('</h>',''))
   h3  = h
   h2  = 10*h
   h1  = 100*h

   xmlTag = dom.getElementsByTagName('rmax')[0].toxml()
   rmax= float(xmlTag.replace('<rmax>','').replace('</rmax>',''))

   Q   = sqrt(d/(d-2))*rs**(d-1)
   M   = r0**d + Q**2/(r0**(d-2))

   v   = sqrt(Rm2 + d**2/4.0)

   return(rs,r0,R2,kx,kt,Rm2,d,h3,h2,h1,rmax,Q,M,v,muq)
