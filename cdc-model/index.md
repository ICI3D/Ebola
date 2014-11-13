---
layout: page
title: CDC's EbolaResponse modeling tool
welcome: An interactive, web-based reimplementation by CAB Pearson
summary: We have re-implemented the model developed by Meltzer et al. This implementation of the model is not endorsed or reviewed by the CDC. Links to source code and license information are provided at the bottom of this page (<a href="#source">link</a>).
---

## Re-Implementation

<iframe name="Model" src="https://pearsonca.shinyapps.io/CDC-Ebola-model/cdc-model.Rmd" width="1100px" height="2200px"></iframe>

The CDC's original implementation of the EbolaResponse modeling tool was published as

> Meltzer, MI _et al_. (2014) Estimating the Future Number of Cases in the Ebola Epidemic --- Liberia and Sierra Leone, 2014-2015. _Morbidity and Mortality Weekly Report_ 63: 1–14. ([Link to paper](http://www.cdc.gov/mmwr/preview/mmwrhtml/su6303a1.htm), [EbolaResponse modeling tool spreadsheet](http://dx.doi.org/10.15620/cdc.24900))

{: #source}
## Source code

This web-based version of their model is written in R and is made available under a Creative Commons attribution (CC-BY) license ([details below](#license)). The code was written by [Dr. Carl A.B. Pearson](http://www.pulliamlab.org/people/pearson), with plots by [Dr. Juliet R.C. Pulliam](http://www.pulliamlab.org/people/pulliam). The source code is available via a public GitHub repository:

> [https://github.com/ICI3D/CDC-Ebola-model](https://github.com/ICI3D/CDC-Ebola-model).

{: #citation}
## Citation

Please cite this reimplementation of the model as:

> <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Pearson, CAB and JRC Pulliam</span>. (2014) Web-based reimplementation of CDC's EbolaResponse modeling tool. [{{ site.subdomainurl }}/cdc-model]({{ site.subdomainurl }}/cdc-model). Source code available at [https://github.com/ICI3D/CDC-Ebola-model](https://github.com/ICI3D/CDC-Ebola-model).

{: #license}
## License

The code is made available under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>. You are free to reuse this code provided that you *give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.* Giving appropriate credit includes citation of this website *and* providing a link to the code repository:

<a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/ICI3D/CDC-Ebola-model" rel="dct:source">https://github.com/ICI3D/CDC-Ebola-model</a>

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />
