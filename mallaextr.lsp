;;;
;;; Andrés Testi 2001.
;;; 
;;; Licensed under the Apache License, Version 2.0 (the "License"); you may not
;;; use this file except in compliance with the License. You may obtain a copy of
;;; the License at
;;; 
;;; http://www.apache.org/licenses/LICENSE-2.0
;;; 
;;; Unless required by applicable law or agreed to in writing, software
;;; distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
;;; WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
;;; License for the specific language governing permissions and limitations under
;;; the License.
;;;


;;;
;;; Extrusor de Mallas en AutoLISP 2000, tomando una spline como guía y otra spline como patrón.
;;; Se invoca por medio del comando MallaExtr.
;;;

;;; Función Producto escalar x Matriz
(defun exA (Esc Matriz)
  (mapcar '(lambda (Fila) (mapcar '(lambda (Elem) (* Esc Elem)) Fila))
	  Matriz
  )
)

;;; Función Matriz Adjunta
(defun MatAdj (Matriz / MatOut FilaOut m n)
  (setq m -1)
  (repeat (length Matriz)
    (setq m (1+ m)
	  n 0
	  FilaOut nil
    )
    (repeat (length (car Matriz))
      (setq FilaOut (append
		      FilaOut
		      (list (* (expt -1 (+ (1+ m) (1+ n)))
			       (Det (mapcar '(lambda (Fila) (Excluir n Fila))
					    (Excluir m Matriz)
				    )
			       )
			    )
		      )
		    )
	    n	    (1+ n)
      )
    )
    (setq MatOut (append MatOut (list FilaOut))
    )
  )
  (Transp MatOut)
)

;;; Función para excluir un atomo de una lista
(defun Excluir (n Matriz / MatOut m)
  (setq m -1)
  (repeat (length Matriz)
    (setq m (1+ m))
    (if	(/= m n)
      (setq MatOut (append MatOut (list (nth m Matriz))))
      MatOut
    )
  )
)


;;; Función Determinate de 3x3
(defun Det3x3 (Matriz / a11 a12 a13 a21 a22 a23 a31 a32 a33)
  (setq	a11 (nth 0 (nth 0 Matriz))
	a12 (nth 1 (nth 0 Matriz))
	a13 (nth 2 (nth 0 Matriz))
	a21 (nth 0 (nth 1 Matriz))
	a22 (nth 1 (nth 1 Matriz))
	a23 (nth 2 (nth 1 Matriz))
	a31 (nth 0 (nth 2 Matriz))
	a32 (nth 1 (nth 2 Matriz))
	a33 (nth 2 (nth 2 Matriz))
  )
  (- (+ (* a11 a22 a33) (* a12 a23 a31) (* a13 a21 a32))
     (+ (* a11 a23 a32) (* a12 a21 a33) (* a13 a22 a31))
  )
)

;;; Función Determinate
(defun Det (Matriz / Determ Adjunto m Fila)
  (setq	m 0
	Determ 0
  )
  (if (= (length Matriz) 1)
    (car (car Matriz))
    (progn
      (foreach Fila Matriz
	(setq Adjunto (mapcar 'cdr (Excluir m Matriz))
	      Elem    (car Fila)
	      m	      (1+ m)
	      Determ
		      (+ Determ (* (expt -1 (+ m 1)) Elem (Det Adjunto)))
	)
      )
    )
  )
)

;;; Función Matriz Inversa
(defun MatInv (Matriz)
  (setq Determ (Det Matriz))
  (if (/= Determ 0)
    (exA (/ 1.0 Determ) (MatAdj Matriz))
  )
)

;;; Función Producto Escalar
(defun ProdEsc (Esc Uvec)
  (mapcar '(lambda (x) (* x Esc)) Uvec)
)

;;; Función Signo
(defun Signo (Numero)
  (if (/= Numero 0)
    (/ (abs Numero) Numero)
    0
  )
)

;;; Función Tangente
(defun Tan (Angulo)
  (/ (Sin Angulo) (Cos Angulo))
)

;;; Función Versor
(defun Versor (Uvec / Unorm)
  (setq	Unorm (Norma Uvec)
  )
  (if (/= Unorm 0)
    (mapcar '(lambda (x) (/ x Unorm)) Uvec)
    Uvec
  )
)

;;; Esta función Suma 2 Vectores Uvec y Vvec
(defun U+V (Uvec Vvec)
  (mapcar '+ Uvec Vvec)
)

;;; Esta función Resta 2 Vectores  Uvec y Vvec
(defun U-V (Uvec Vvec)
  (mapcar '- Uvec Vvec)
)

;;; Esta función calcula la norma de un Vector Uvecn
(defun Norma (Uvec)
  (sqrt (apply '+ (mapcar '(lambda (x) (expt x 2.0)) Uvec)))
)

;;; Esta función calcula el producto punto de los Vectores Uvecp y Vvecp
(defun U*V (Uvec Vvec)
  (apply '+ (mapcar '* Uvec Vvec))
)

;;; Esta función calcula el producto cruz entre 2 vectores de R3
(defun UxV (Uvec Vvec / U1 U2 U3 V1 V2 V3)
  (setq	U1 (car Uvec)
	U2 (cadr Uvec)
	U3 (caddr Uvec)
	V1 (car Vvec)
	V2 (cadr Vvec)
	V3 (caddr Vvec)
  )
  (list	(- (* U2 V3) (* U3 V2))
	(- (* U3 V1) (* U1 V3))
	(- (* U1 V2) (* U2 V1))
  )
)

;;; Esta función multiplica una Matriz por un Vector
;;; Los Argumentos Mat y Vec son las listas que contienen a la Matriz y al
;;; Vector respectivamente
(defun AxU
	   (Mat Vec)
  (mapcar '(lambda (Matfila) (apply '+ (mapcar '* Matfila Vec)))
	  Mat
  )
)

;;; Función Matriz Transpuesta
(defun Transp (Matriz / Mresto Mnueva)
  (setq Mresto Matriz)
  (repeat (length (car Matriz))
    (setq
      Mnueva (append Mnueva (list (mapcar 'car Mresto)))
      Mresto (mapcar 'cdr Mresto)
    )
  )
  Mnueva
)


;;; Función Producto Matricial
(defun AxB (MatrizA MatrizB / FilaA ColB RestoB ElemC FilaC MatrizC)
  (foreach FilaA MatrizA
    (setq RestoB MatrizB
	  FilaC	 nil
    )
    (repeat (length (car MatrizB))
      (setq ColB   (mapcar 'car RestoB)
	    RestoB (mapcar 'cdr RestoB)
	    ElemC  (apply '+ (mapcar '* FilaA ColB))
	    FilaC  (append FilaC (list ElemC))
      )
    )
    (setq MatrizC (append MatrizC (list FilaC)))
  )
)

;;; Esta función calcula el coseno del ángulo interno entre 2 vectores
;;; Los argumentos Uvecc y Vvecc son las listas que contienen las  coordenadas
;;; de los Vectores U y V respectivamente.
(defun CosUV (Uvec Vvec / Unorm Vnorm)
  (setq	Unorm (Norma Uvec)
	Vnorm (Norma Vvec)
  )
  (if (and (/= Unorm 0) (/= Vnorm 0))
    (/ (U*V Uvec Vvec) (* Unorm Vnorm))
    1.0
  )
)

;;; Esta función calcula el seno del ángulo interno entre 2 vectores de R2
;;; Los argumentos Uvecc y Vvecc son las listas que contienen las  coordenadas
;;; de los Vectores U y V respectivamente.
(defun senUV (Uvecs Vvecs)
  (sqrt	(- 1.0 (expt (cosuv uvecs vvecs) 2.0))
  )
)

;;; Esta función calcula el seno del ángulo interno entre 2 vectores de R2
;;; Los argumentos Uvecc y Vvecc son las listas que contienen las  coordenadas
;;; de los Vectores U y V respectivamente.
;;; El signo del seno dependerá del vector que se halla ingresado primero
(defun SenUV (Uvecs Vvecs / Unorms Vnorms Uvecs1 Uvecs2 Vvecs1 Vvecs2)
  (setq	Unorms (Norma Uvecs)
	Vnorms (Norma Vvecs)
	Uvecs1 (car Uvecs)
	Uvecs2 (cadr Uvecs)
	Vvecs1 (car Vvecs)
	Vvecs2 (cadr Vvecs)
  )
  (cond	((or (= Unorms 0) (= Vnorms 0)) 0.0)
	((= Vvecs2 0) (/ Uvecs2 Unorms))
	(t
	 (/ (- (* (/ Vnorms Unorms) Uvecs1)
	       (* (cosuv Uvecs Vvecs) Vvecs1)
	    )
	    Vvecs2
	 )
	)
  )
)

;;; Esta función calcula el seno del ángulo interno entre 2 vectores de R2
;;; Los argumentos Uvecc y Vvecc son las listas que contienen las  coordenadas
;;; de los Vectores U y V respectivamente.
;;; El signo del seno dependerá del vector que se halla ingresado primero
(defun SenUV2 (Uvecs Vvecs / Unorms Vnorms Uvecs1 Uvecs2 Vvecs1 Vvecs2)
  (setq	Unorms (Norma Uvecs)
	Vnorms (Norma Vvecs)
	Uvecs1 (car Uvecs)
	Uvecs2 (cadr Uvecs)
	Vvecs1 (car Vvecs)
	Vvecs2 (cadr Vvecs)
  )
  (cond	((or (= Unorms 0) (= Vnorms 0)) 0.0)
	((= Uvecs2 0) (/ Vvecs2 Vnorms))
	(t
	 (/ (- (* (/ Unorms Vnorms) Vvecs1)
	       (* (cosuv Uvecs Vvecs) Uvecs1)
	    )
	    Uvecs2
	 )
	)
  )
)

;;; Esta  función  transforma por Traslación una  lista de puntos (dados en
;;; vectores  de  R3)  ListaInp, según  dos vectores de posición, InpVec de
;;; entrada, y OutVec de salida
(defun TransTrasLista (Uvec Vvec ListaPuntos)
  (mapcar '(lambda (Punto) (mapcar '+ (mapcar '- Vvec Uvec) Punto))
	  ListaPuntos
  )
)

;;; Esta función transforma por Traslación un Vector posición en R3, Wvec, según
;;; dos vectores de posición, InpVec de entrada, y OutVec de salida
(defun TransTrasVec (InpVec OutVec Wvec)
  (U+V Wvec
       (U-V OutVec InpVec)
  )
)

;;; Esta Función Activa el VLisp
(vl-load-com)


;;; Esta función permite que el Usuario seleccione las curvas y devuelve una lista
;;; conteniendo la denominación de ambas curvas
(defun SelCurvas (/ Patron Pista)
  (while (not
	   (setq Patron (entsel "\nSeleccione el patrón a extruir:"))
	 )
  )
  (while (not (setq Pista (entsel "\nSeleccione la Pista de extrusión:")
	      )
	 )
  )
  (list	(car Patron)
	(car Pista)
	(cadr Pista)
  )
)

;;; Esta función devuelve la lista de coordenadas de los puntos que resultan de dividir una
;;; curva Curvf en una cantidad determinada de segmentos Segmf
(defun Fracmentar
		  (Curvf      Segmf	 /	    Buclef
		   Listaf     ComCurva	 FinCurva   Distancia
		   LongSeg
		  )
  (setq	ComCurva  (vlax-curve-getStartParam Curvf)
	FinCurva  (vlax-curve-getendParam Curvf)
	LongSeg	  (/ (vlax-curve-getdistatParam Curvf FinCurva) Segmf)
	Buclef	  0
	Distancia 0
	Listaf	  nil
  )
  (while (< Buclef Segmf)
    (setq Listaf
		    (append Listaf
			    (list (vlax-curve-getPointAtdist Curvf Distancia))
		    )
	  Buclef    (1+ Buclef)
	  Distancia (* Buclef LongSeg)
    )
  )
  (append Listaf
	  (list (vlax-curve-getPointAtParam Curvf FinCurva))
  )
)

;;; Esta es la función principal de malla extruida
(defun C:MallaExtr (/		     ListaCurvas      Patron
		    Pista	     PuntoOrigen      Tab1
		    Tab2	     OrigenPista      ListaPatronWCS
		    ListaPistaWCS    DerivInicial     BucleME
		    ListaPatronOCS   ListaPistaOCS    ListaMallaWCS
		    ListaMallaOCS    PuntoPistaOCS    PuntoPistaWCS
		   )
  (setq	ListaCurvas (Selcurvas)
	Patron	    (car Listacurvas)
	Pista	    (cadr ListaCurvas)
	PuntoOrigen (caddr ListaCurvas)
	Tab1	    (- (getvar "surftab1") 1)
	Tab2	    (- (getvar "surftab2") 1)
  )
  (vlax-ename->vla-object Patron)
  (vlax-ename->vla-object Pista)
  (setq	OrigenPista
	 (vlax-curve-getStartPoint Pista)
	ListaPatronWCS
	 (Fracmentar Patron Tab1)
	ListaPistaWCS
	 (Fracmentar Pista Tab2)
	DerivInicial
	 (vlax-curve-getFirstDeriv
	   Pista
	   (vlax-curve-getStartParam Pista)
	 )
	BucleME	0
	ListaPatronOCS
	 (TransTrasLista OrigenPista '(0 0 0) ListaPatronWCS)
	ListaPistaOCS
	 (TransTrasLista OrigenPista '(0 0 0) ListaPistaWCS)
	ListaMallaWCS
	 nil
	ListaMallaOCS
	 nil
  )

  (foreach PuntoPistaOCS ListaPistaOCS
    (setq PuntoPistaWCS	(nth BucleME ListaPistaWCS)
	  ListaMallaOCS
			(append	ListaMallaOCS
				(TransTrasLista
				  '(0 0 0)
				  PuntoPistaOCS
				  (TransRotLista
				    DerivInicial
				    (vlax-curve-getFirstDeriv
				      Pista
				      (vlax-curve-getparamatpoint Pista PuntoPistaWCS)
				    )
				    ListaPatronOCS
				  )
				)
			)
	  BucleME
			(1+ BucleME)
    )
  )
  (setq	ListaMallaWCS
	 (TransTrasLista '(0 0 0) OrigenPista ListaMallaOCS)
  )
  (TrazarMallaDXF Tab2 Tab1 ListaMallaWCS)
)

;;; Esta función traza una malla DXF, según m y n y una lista de puntos
(defun TrazarMallaDXF
		      (Mseg Nseg ListaMalla / Vertice)
  (entmake (list (cons 0 "POLYLINE")
		 (cons 70 16)
		 (cons 71 (1+ Mseg))
		 (cons 72 (1+ Nseg))
		 (cons 66 1)
	   )
  )
  (foreach Vertice ListaMalla
    (entmake (list (cons 0 "VERTEX")
		   (cons 100 "AcDbVertex")
		   (cons 10 Vertice)
		   (cons 70 64)
	     )
    )
  )
  (entmake (list (cons 0 "SEQEND")))
)


(defun Prueba ()
  (fracmentar (vlax-ename->vla-object (car (entsel))) 7)
)

(defun Prueba2 (/ ListaPuntos Curva Punto Derivada Vector)
  (setq	Curva	    (vlax-ename->vla-object (car (entsel)))
	ListaPuntos
		    (fracmentar Curva 10)
  )
  (foreach Punto ListaPuntos
    (setq Derivada   (vlax-curve-getFirstDeriv
		       Curva
		       (vlax-curve-getparamatpoint Curva Punto)
		     )
	  Vector     (Prodesc 2 (Versor Derivada))
	  PuntoFinal (U+V Punto Vector)
    )
    (entmake (list (cons 0 "CIRCLE")
		   (append '(10) Punto)
		   (cons 40 0.1)
	     )
    )
    (entmake (list (cons 0 "LINE")
		   (cons 100 "AcDbLine")
		   (append '(10) Punto)
		   (append '(11) PuntoFinal)
	     )
    )
  )
)


;;; Transformación de Rotación
(defun TransRotLista (Uvec     Vvec	ListaPuntos	  /
		      ir       jr	kr	 vr	  SenRot
		      CosRot   MatRot	MatTrans MatTInv  MaTFin
		     )
  (if (not (equal (versor Uvec)(versor Vvec)))
     (progn
       (setq ir	      (Versor Uvec)
	     kr	      (Versor (UxV Uvec Vvec))
	     jr	      (Versor (UxV kr Uvec))
	     MatTrans (Transp (list ir jr kr))
	     MatTInv  (MatInv MatTrans)
	     vr	      (Versor (AxU MatTInv Vvec))
	     CosRot (car vr)
	     SenRot (cadr vr)
	     MatRot   (list (list CosRot (* -1 SenRot) 0)(list SenRot CosRot 0)'(0 0 1))
	     MatFin   (AxB MatTrans (AxB MatRot MatTInv))
       )
       (mapcar '(lambda (Punto) (AxU MatFin Punto)) ListaPuntos)
     )
     ListaPuntos
  )
)
