# Capteur MEMS résonant pour la mesure de masse ultra-précise

Ce projet modélise et simule le comportement d’un capteur MEMS résonant pour la mesure de très faibles masses, à l’aide de méthodes numériques avancées.

## Description

Le capteur étudié est basé sur une nano-poutre soumise à une force électrostatique. L’ajout d’une masse modifie la fréquence de résonance du système, ce qui permet de quantifier la masse déposée. Deux méthodes numériques sont explorées :

- **Mesure du décalage en fréquence** : suivi de la variation de la fréquence de résonance en fonction de la masse ajoutée.
- **Détection par saut d’amplitude** : détection d’un saut d’amplitude lors d’un balayage en fréquence, permettant de détecter des masses encore plus faibles.

Les résultats montrent la faisabilité de ces approches et ouvrent la voie à la conception de capteurs MEMS/NEMS de très haute sensibilité.

## Structure du projet

- `newmark.py`, `newmark_masse_ajoutee.py`, `newmark_saut.py` : modules de calcul et d’intégration numérique (méthode de Newmark, courbes de réponse, etc.)
- `trace_delta_m_f.py` : analyse du décalage en fréquence en fonction de la masse ajoutée
- `saut_final.py` : simulation du saut d’amplitude
- `Courbes de réponse/` : données de courbes de réponse pré-calculées
- `jupyter-MEMS.ipynb` : rapport interactif et visualisations

## Installation

1. Cloner le dépôt
2. Installer les dépendances Python :

   ```bash
   pip install -r requirements.txt
   ```

## Utilisation

- Lancer les scripts Python pour générer les courbes de réponse et les analyses.
- Ouvrir le notebook `jupyter-MEMS.ipynb` pour explorer le rapport interactif et les visualisations.

## Dépendances principales

Voir `requirements.txt` pour la liste complète.

## Auteurs

- Thomas Bertrand
- Marc Doumit
- Taihani Chan
- Tristan Hausermann

Référent : Sébastien Baguet
