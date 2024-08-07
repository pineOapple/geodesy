## 2.1 Bewegung der Sonne

Überlegen Sie selber die Rektaszension $\alpha$ und Deklination $\delta$ der Sonne am Sommer-, Winter- und Frühlingsanfang

1. Die Definition des Frühlingspunktes lautet: Nullstelle der Deklination mit positiver Steigung, also gerade die Position der Sonne am 21.03. um UT1=12h bei $\Lambda,\Phi=0°$ unter Vernachlässigung der Präzession. 
2. Der Frühlingspunkt und die Knotenlinie der Erde sind zu diesem Zeitpunkt kollinear. Und wir erhalten $\alpha, \delta=0$ für den **Frühlingsanfang am 21.03**
3. Mit der Umlaufzeit der Erde von etwa 365 Tagen und der Annahme einer kreisförmigen Umlaufbahn erwarten wir nach einem halben Jahr eine Rektaszension von $\alpha=180°$ und wieder sind Frühlingspunkt und Knotenlinie kollinear, also $\delta=0°$ zum **Herbstanfang am 21.09.**.
4. Die Sonne steht im unseren Breitengrad von ca. 48° auf der Nordhalbkugel und unserer Nähe am Nullmeridian zum **Sommeranfang am 21.06.** am höchsten, also eine positive Deklination von $\delta=23.5°$ und einer aus der allgemeinen Definition der Richtung der Umlaufbahn der Erde eine Rektaszension $\alpha=90°$.
5. Übrig bleibt der **Winteranfang am 21.12.** mit einer negativen Deklination von $\delta=-23.5°$ und $\alpha=270°$

und wandeln Sie $\alpha$ und $\delta$ in einen kartesischen Vektor $r_i'$ um. 

1. Aus der Definition der sphärischen Koordinaten an der Einheitskugel finden wir:
$$
r_i' = \begin{pmatrix}
\cos \delta \cos \alpha \\
\cos \delta \sin \alpha \\
\sin \delta
\end{pmatrix}
$$

Die Koordinaten $\Phi$ und $\Lambda$ dürfen Sie selbst für Ihren Lieblingsort auswählen. 

1. Aus Gründen der Verifizierbarkeit unserer Simulationen wählen wir zunächst die Koordinaten eines Atoms einer Boje bei $\Phi=0.000000000000000°$ und $\Lambda=0.000000000000000°$. Die geodätischen und geozentrischen Koordinaten unterscheiden sich an diesem Punkt nicht. 
2. Wir führen noch die Koordinaten der Neumayer III Station bei $\Phi=-70.6374°$ und $\Lambda=8.2609°$ auf der Antarktis zum Vergleich ein.

Ihre Aufgabe ist es, die Position der Sonne im lokalen System zu berechnen, also Azimut $A$ und Zenitwinkel $z$. Die Transformation lautet vollständig:
$$ r_g' = S_1 (R_z(90^\circ - \phi)) R_3 (\alpha_{\text{GAST}}) N P r_i' $$
mit
$$ N = R_1 (-\epsilon - \Delta \epsilon) R_3 (-\Delta \psi) R_1 (\epsilon) $$
$$ P = R_3 (-z) R_2 (\theta) R_3 (-\xi) . $$
1. Aufgrund der Definition der internen Drehmatrizen in *matlab* kehren wir alle Vorzeichen der Argumente in Drehmatrizen um und verwenden `rotx`, `roty` und `rotz`.
2. Sie finden den zugehörigen Code im Anhang `skyboxxer.m`

Zeichnen Sie mittels `skyplot` den scheinbaren Umlauf der Sonne während eines ganzen Tages. Wie lang ist der Tag (Sonnenaufgang bis -untergang) an diesen beiden Tagen?
![[Pasted image 20240720135917.png]]

![[Pasted image 20240720140015.png]]
![[Pasted image 20240720140005.png]]

![[Pasted image 20240720140039.png]]


```
Kein Sonnenuntergang für: 21/January/2024 verzeichnet

Sonnenaufgang   für 21/February/2024: 4h 42min
Sonnenuntergang für 21/February/2024: 20h 48min

Sonnenaufgang   für 21/March/2024: 6h 54min
Sonnenuntergang für 21/March/2024: 18h 30min

Sonnenaufgang   für 21/April/2024: 9h 6min
Sonnenuntergang für 21/April/2024: 15h 54min

Kein Sonnenaufuntergang für: 21/May/2024 verzeichnet
Kein Sonnenaufuntergang für: 21/June/2024 verzeichnet
Kein Sonnenaufuntergang für: 21/July/2024 verzeichnet

Sonnenaufgang   für  21/August/2024: 9h 6min
Sonnenuntergang für  21/August/2024: 16h 12min

Sonnenaufgang   für  21/September/2024: 6h 36min
Sonnenuntergang für  21/September/2024: 18h 24min

Sonnenaufgang   für  21/October/2024: 4h 6min
Sonnenuntergang für  21/October/2024: 20h 36min

Kein Sonnenuntergang für: 21/November/2024 verzeichnet
Kein Sonnenuntergang für: 21/December/2024 verzeichnet
```

Diskussion.

1. Mithilfe der Inversen unserer Transformation können wir das Koordinatensystem wieder auf ein "klassisches" Koordinatensystem legen. Damit waren die 3D plots möglich und auch die Verifikation für korrekte Transformation und dem korrekten Umgang mit den Daten.
2. Für das **Frühlingsäquinoktium** erhalten wir eine **Tageslänge** von 11h und 36m
3. Für die **Sommersonnenwende** erhalten wir eine **Tageslänge** von 0h, die Sonne geht hier an der Antarktis nicht auf
4. Für das **Herbstäquinoktium** erhalten wir eine **Tageslänge** von 11h und 46m
5. Für die Wintersonnenwende erhalten wir eine Tageslänge von 
```octave
% [4] Create inertial transformation matrix and rotate GWC
while(true)
	for j = 1:0.1:24
		[R, G] = GDS_INT_TO_LCL(lambda,phi,yr,m,d,j);
		gw_c = G*gw_i;
		set(h, 'XData', gw_c(1, :), 'YData', gw_c(2, :), 'ZData', gw_c(3, :));
		pause(0.05);
	end
end
```
Zur Berechnung von Nutationswinkeln $\{\epsilon, \Delta \epsilon, \Delta \psi\}$ und Präzessionswinkeln $\{\zeta, \theta, \xi\}$ werden den Ihnen die m-files nutwink und prezwink zur Verfügung gestellt. Zusätzlich werden hier julianjh und jul2gast benötigt.
