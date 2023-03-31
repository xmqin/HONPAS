project: SIESTA
version: {!../../version.info!} 
author: SIESTA Group
src_dir: ../../Src
         ../../Util/Bands
         ../../Util/COOP
output_dir: ./doc-ford-build
page_dir: ./ford-pages
media_dir: ./ford-media
display: public
         private
source: true
search: false
graph: false
exclude_dir: ../../Src/MPI
             ../../Src/fdf
             ../../Src/fdict
             ../../Src/ncdf
             ../../Src/Libs
             ../../Src/Orphans
             ../../Src/SiestaXC
             ../../Src/wxml
             ../../Src/xmlparser
exclude:   m_gauss_fermi_17.f90
		   m_gauss_fermi_18.f90
		   m_gauss_fermi_19.f90
		   m_gauss_fermi_20.f90
		   m_gauss_fermi_22.f90
		   m_gauss_fermi_24.f90
		   m_gauss_fermi_26.f90
		   m_gauss_fermi_28.f90
		   m_gauss_fermi_30.f90
		   m_gauss_fermi_inf.f90
		   spinorbit.f
project_website: https://www.icmab.es/siesta
summary: ![SIESTA](logo) <br/>
         A first-principles materials simulation code 
         using Density Functional Theory.
         {: style="text-align: center" }
preprocess: true
preprocessor: gfortran -E -P
docmark_alt: *
predocmark: >
fpp_extensions: F90
                F
license: gfdl
extra_filetypes: sh #
                 inc ! 
md_extensions: markdown.extensions.toc

@note
This is an early stage work-in-progress build of developers documentation for SIESTA.
@endnote

## Project Dashboard

@todo
Add more topics
@endtodo


<!-- useful options -->
<!-------------------->
<!-- graph_maxdepth: 100 -->
<!-- graph_maxnodes: 100 -->


<!-- predocmark: >  -->
