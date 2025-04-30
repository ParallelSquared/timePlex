
![Logo](timePlex_figure.png)

## Article: [Increasing mass spectrometry throughput using time-encoded sample multiplexing (Derks et al, 2025)](https://www.biorxiv.org/)

 

<h2 style="letter-spacing: 2px; font-size: 26px;" id="RAW-data">

Demonstrated aim:

</h2>

This work demonstrates an approach to multiplex samples for proteomic analysis using temporal encoding, which we refer to as 'timePlex'. This entailed the development of experimental implementations and computational analyses to deconvolve time-encoded samples, as integrated as a module in [JMod](https://github.com/ParallelSquared/jmod). Experiments were performed on 8-12 ng standards to benchmark proteomic coverage and quantitative accuracy of this new 'timePlex' data type to established methods (LF non-multiplexed DIA and plexDIA). Single-cell data was also acquired using combinatorial 3x3 timePlex and plexDIA, enabling a 9-plex.

 

<h2 style="letter-spacing: 2px; font-size: 26px;" id="RAW-data">

Methods:

</h2>

#### Bulk data:

Proteomics methods: [plexDIA](https://scp.slavovlab.net/plexDIA) & [timePlex](https://www.parallelsq.org/technology-and-science/timePlex)<br>

 

#### Single-cell data:

Proteomics methods: [plexDIA](https://scp.slavovlab.net/plexDIA) & [timePlex](https://www.parallelsq.org/technology-and-science/timePlex)<br>

Sample preparation method: [nPOP](https://scp.slavovlab.net/nPOP)<br>  

<h2 style="letter-spacing: 2px; font-size: 26px;" id="plexDIA-data">

Data:

</h2>

All raw and processed data from the [article](https://www.biorxiv.org/) are organized in this MassIVE repository: [MSV000097736](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=7193ea0d007741c680f22ec005718e2b).

 

<h2 style="letter-spacing: 2px; font-size: 26px;" id="code">

Code:

</h2>

[JMod](https://github.com/ParallelSquared/jmod) was used for searching all data. The exact version of JMod which was used for searching these raw data is available for download [here](https://github.com/ParallelSquared/timePlex/blob/main/JMod_timePlex.zip).

The [{targets} R package](https://books.ropensci.org/targets/) was used to ensure repeatability of our downstream analyses. R scripts corresponding to benchmarking and single-cell analyses can be found in the ["Code" folder of our Github](https://github.com/SlavovLab/Single-nucleus-proteomics/tree/main/Code). The raw data, processed data, libraries, and meta data required to repeat these analyses can be found at our MassIVE repository: [MSV000097736](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=7193ea0d007741c680f22ec005718e2b)

 

 

 
