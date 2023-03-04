
# DAFE
## Description
Cette application Rshiny est composé en 2 partie : Analyse différentielle et Enrichissement fonctionnelle. 
### Analyse differentielle
La première partie prend en entrée des fichiers de comptage brut et réalise la normalisation et l'analyse différentielle de ces derniers.
Les fichiers de résultats sont les comptages normalisés et le fichier d'analyses différentielle. De differentes représenattion graphique sont également disponible( Heatmpa des réplicas, Volcano plot ...)
### Enrichissement fonctionnel
Cette secone partie prend en entier un fichier d'analyse différentielle, et va estimer par deux méthodes (ORA et GSEA) si il existe un enrichissement de certains GO terms, pathway KEGG ou Reactome. Pour chaque enrichissement un fichier de sortie est disponible, composé entre autre de la stat de test, de la pvalue, et pvalue ajustée. De plus plusieurs représenattion grapique sont également disponible : Barplot, Goplot, Dotplot, GSEAplot ...


## Auteur

- [@Audrey](https://www.github.com/audrey2)


## Installation



#### Installer DAFE

```bash
  git clone https://github.com/audrey2/Differential_Analaysis_Rshiny.git

```
#### Lancer DAFE 
```bash
  cd Differential_Analaysis_Rshiny
  Rscript -e 'golem::run_dev()'
```

#### Installer les packages necessaires à DAFE
```bash
  Rscript install_packages_necessary.R
```


## Demonstration

Il existe des données test inclus directement dans l'application, ainsi une demo de l'application est disponible en cliquant sur le bouton 'use data test' du second onglet:


![image](https://github.com/audrey2/DAFE/blob/main/Data_test_DAFE.png?raw=true)


