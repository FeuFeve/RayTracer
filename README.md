# Programmation 3D - Projet de Raytracing


#### Auteur : Clément POTIN, M1 IMAGINA


J'ai codé sous l'IDE CLion sous Windows, et n'ai donc pas eu à compiler manuellement mon projet.

Des images de rendus obtenus lors de ma progression sont disponibles dans le dossier
`cmake-build-debug/Progress`, et de bugs que j'ai trouvé intéressant à conserver dans le dossier
`cmake-build-debug/Bugs`.


## Fonctionnement

Je n'ai conservé que la scène de la boîte de Cornell sur les 3 scènes d'origine, et j'ai ajouté une
boîte de Cornell dans laquelle les sphères sont remplacées par des objets composés de triangles.
Appuyer sur la touche `1` permet donc d'accéder à la scène d'origine, et sur `2` permet d'accéder à
la scène au style "low poly".

Pour effectuer un rendu, j'ai laissé la touche `R`.
Il est également possible de cliquer sur la fenêtre pour lancer un CastRayDebug, que j'ai modifié
pour obtenir certaines informations qui m'ont été très utile pour débugger tout au long du projet.


## Détails du code

Les rendus sont effectués via un modèle de Phong (lumières ambiante + diffuse + spéculaire (qui est
ensuite devenue la refléxion/les miroirs)).
N'ayant pas trouvé d'informations précises à ce sujet, j'ai pris la liberté de :

	- Définir comme couleur minimum la couleur de l'objet courant * son facteur de réflexion
	ambiante. Un objet RGB = (200, 100, 40) dont le coefficient Ka est à 0.1 ne pourra donc pas
	être plus sombre que RGB = (20, 10, 4), même si le pixel est entièrement dans l'ombre ;

	- Définir la priorité de la lumière sur la réflexion, puis sur la réfraction. Autrement dit :
	`color = Ks * reflectedColor + (1-Ks) * (Kt * refractedColor + (1-Kt) * color)`
	Si Ks = 1 (miroir parfait), alors il n'y aura aucune réfraction, même si la transparence (Kt)
	est à 1 (verre parfait).
	De cette façon je peux obtenir des effets assez réalistes sur des sphères en verre légèrement
	teintées, en leur donnant des coefficients tels que Ks = 0.1, Kt = 0.8, Kr = 1.5 (ici la
	couleur proviendra à 10% de la réflexion, à 0.8 * 90% = 72% de la réfraction, et aux 18%
	restants de la couleur d'origine de la sphère).

	- Refaire le système d'ombres d'une façon qui me paraissait plus intuitive.


Pour chercher des objets composés de triangle, j'ai ajouté un type d'objet (`Object::MeshObject`),
qui a pour rôle de charger le modèle (fichier .obj) et de réaliser les calculs d'intersections
rayon/triangle, ainsi que d'intersection rayon/bounding box (une optimisation que j'ai implémenté
pour réduire le nombre de tests rayon/triangle réalisés).

Le calcul d'intersection rayon/triangle est l'implémentation C++ de l'Algorithme d'intersection de
Möller–Trumbore, réadapté pour le projet, et dont la version de base est disponible ici :
https://fr.wikipedia.org/wiki/Algorithme_d%27intersection_de_M%C3%B6ller%E2%80%93Trumbore

Quant à l'optimisation "rayon/bounding box" utilise les vec3 "box_min" et "box_max" de la mesh des
objets. Les fonctions `Object::raySquareIntersection` (et non pas `Square::raySquareIntersection`,
qui utilise forcément un objet "Square") et `MeshObject::rayBoundingBoxIntersection` servent à
cela. Le principe est de vérifier si le rayon entre dans la bouding box de l'objet (en réalisant 6
intersections rayon/rectangle) et de ne réaliser les intersections rayon/triangle de l'objet que si
c'est le cas.

Avec les paramètres suivants :
	- nbPoints = 50 (main.cpp, ligne 37)
	- maxDepth = 5 (main.cpp, ligne 40)
Et sur la seconde scène, le rendu prend sur ma machine 376.15s sans utiliser l'optimisation, contre
seulement 133.25s en l'utilisant (pour tester sans : commenter la ligne 179 de `Object.cpp`
(`if (rayBoundingBoxIntersection(p0, V))`)).


Enfin, j'ai également implémenté une classe `Chrono.cpp`, très utile pour afficher des temps de
rendu ou calculer des ETA (temps estimé avant la fin du rendu). Leur affichage se fait
automatiquement dans la console au moment du rendu.


## Points manquants au projet

La seule partie que je n'ai pas conservée (car incomplète) est le dernier point de la Phase 3 :
"Interpoler les normales et calculer l’ombrage et la couleur à l’aide celle-ci" (je me suis dit que
l'optimisation réalisée compenserait).

Les autres points devraient être complétés et entièrement fonctionnels.


## Liens utiles

J'ai utilisé Git et GitHub pour travailler sur un répertoire privé (que je rendrai public le 26
mai).

Lien : https://github.com/FeuFeve/RayTracer
Images de progression : https://github.com/FeuFeve/RayTracer/tree/master/cmake-build-debug/Progress
Images de différents bugs : https://github.com/FeuFeve/RayTracer/tree/master/cmake-build-debug/Bugs
