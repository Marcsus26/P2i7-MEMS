# Visualisation des bassins d'attraction et analyse du saut d'amplitude dans un système MEMS non linéaire

Ce projet propose des outils de simulation et de visualisation pour l'étude des phénomènes non linéaires dans les capteurs MEMS, en particulier l'analyse des bassins d'attraction et du saut d'amplitude lors de l'ajout de masse.

## Description

Le code permet de :

- Simuler la réponse dynamique d'un résonateur MEMS soumis à une excitation harmonique.
- Visualiser les bassins d'attraction et les courbes de réponse (amplitude/fréquence) pour différentes conditions (avant et après ajout de masse).
- Animer l'évolution du système dans un notebook Jupyter, avec une gestion correcte de l'affichage pour éviter les figures blanches ou clignotantes.

Les principaux fichiers et modules :

- `jupyter-MEMS.ipynb` : Notebook principal pour la visualisation et l'animation.
- `saut_final.py` : Script d'analyse du saut d'amplitude.
- `newmark.py`, `newmark_saut.py`, `newmark_masse_ajoutee.py` : Modules d'intégration numérique (méthode de Newmark) adaptés au système MEMS.
- `trace_delta_m_f.py` : Outils de tracé pour l'analyse de la variation de masse.
- `Courbes de réponse/` : Dossier contenant les courbes de réponse pré-calculées pour différentes valeurs de masse.
  
## Utilisation

1. **Installation des dépendances**

   Installez les dépendances Python nécessaires :

   ```bash
   pip install -r requirements.txt
   ```

## Utilisation

1. **Installation des dépendances**
   Installez les dépendances Python nécessaires :

   ```bash
   pip install -r requirements.txt
   ```

2. **Lancement des simulations**
   - Ouvrez le notebook `jupyter-MEMS.ipynb` pour explorer les animations et visualisations interactives.
   - Utilisez `saut_final.py` pour générer et analyser les courbes de saut d'amplitude.

3. **Organisation des données**
   Les courbes de réponse doivent être placées dans le dossier `Courbes de réponse/` avec la nomenclature attendue (voir le notebook ou le script pour les noms de fichiers).

## Dépendances principales

- `numpy` : Calcul numérique
- `matplotlib` : Visualisation et animation
- `numba` : Accélération des calculs (JIT)
- `ipython` : Affichage interactif dans les notebooks (clear_output, display)

## Licence

Ce projet est distribué sous licence MIT.
