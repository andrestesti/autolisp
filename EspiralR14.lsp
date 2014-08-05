;;;
;;; Andrés Testi 2000.
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
;;; Generador de hélices en AutoLISP R14.
;;; Se invoca por medio del comando ESPiral.
;;;

;;; DEFINICION DE VARIABLES GLOBALES
(setq Espvuelta 1)
(setq Espaltura 0)


;;; FUNCION DE ECUACIONES PARAMETRICAS
(defun ESP:Parametrica
       (Angulo ESP-1 ESP-2 ESP-3 EscH EscV Ainicial Rinicial)
  (entmake
    (list
      (cons 0 "VERTEX")
      (cons 100 "AcDb3dPolylineVertex")
      (cons
	10
	(trans
	  (list
	    (+ ESP-1
	       (* (+ (* EscH (- Angulo Ainicial)) Rinicial) (cos Angulo))
	    )
	    (+ ESP-2
	       (* (+ (* EscH (- Angulo Ainicial)) Rinicial) (sin Angulo))
	    )
	    (+ ESP-3
	       (* EscV (- Angulo Ainicial))
	    )
	  )
	  1
	  0
	)
      )
      (cons 70 32)
    )
  )
)


;;; FUNCION PARA CALCULAR LA ESPIRAL
(defun ESP:Calculo (ListaCalcular     /	       ESP-H	ESP-K
		    ESP-L    ESP-X    ESP-Y    ESP-Z	Radio0
		    Radio1   Beta     Beta0    Beta1	EHori
		    EVert    PI/12    Menor    Mayor
		   )
  (setq	ESP-H  (cdr (assoc 1 ListaCalcular))
	ESP-K  (cdr (assoc 2 ListaCalcular))
	ESP-L  (cdr (assoc 3 ListaCalcular))
	Radio0 (cdr (assoc 10 ListaCalcular))
	Radio1 (cdr (assoc 20 ListaCalcular))
	Beta0  (cdr (assoc 120 ListaCalcular))
  )
  (setq	EHori (/ (- Radio1 Radio0) (* Espvuelta (* 2 pi)))
	EVert (/ Espaltura (* Espvuelta (* 2 pi)))
	PI/12 (/ pi 12)
	Beta1 (+ Beta0 (* 2 (* pi Espvuelta)))
	Menor (min Beta0 Beta1)
	Mayor (max Beta0 Beta1)
	Beta  (if (< Menor 0)
		(* (fix (/ Menor PI/12)) PI/12)
		(* (+ (fix (/ Menor PI/12)) 1) PI/12)
	      )
  )
  (entmake (list (cons 0 "POLYLINE")
		 (cons 70 8)
		 (cons 75 5)
		 (cons 66 1)
	   )
  )
  (ESP:Parametrica
    Menor ESP-H	ESP-K ESP-L EHori EVert	Beta0 Radio0)
  (while (< Beta Mayor)
    (ESP:Parametrica
      Beta ESP-H ESP-K ESP-L EHori EVert Beta0 Radio0)
    (setq Beta (+ Beta PI/12))
  )
  (ESP:Parametrica
    Mayor ESP-H	ESP-K ESP-L EHori EVert	Beta0 Radio0)
  (entmake (list (cons 0 "SEQEND")))
)


;;; FUNCION TOMA DE PUNTOS
(defun ESP:TomaPuntos
		      (/ Centro Radio0 Radio1 Altura Vueltas Ceros)
  (initget 1)
  (setq Centro (getpoint "\nCentro de la Espiral:"))
  (if
    (setq
      Radio0
       (getpoint Centro (strcat "\nPunto inicial <centro>:"))
    )
     t
     (setq Radio0 Centro)
  )
  (if (setq Radio1 (getdist Centro "\nRadio 2 <hélice>:"))
    t
    (setq
      Radio1 (distance
	       Radio0
	       (list (car Centro) (cadr Centro) (caddr Radio0))
	     )
    )
  )
  (setq Ceros (getvar "DIMZIN"))
  (setvar "DIMZIN" 1)
  (if (setq Altura (getdist
		     (strcat "\nAltura <" (rtos Espaltura) ">:")
		   )
      )
    (setq Espaltura Altura)
  )
  (initget 2)
  (if (setq Vueltas (getreal (strcat "\nNúmero de vueltas <"
				     (rtos Espvuelta)
				     ">:"
			     )
		    )
      )
    (setq Espvuelta Vueltas)
  )
  (setvar "DIMZIN" Ceros)

  (list
    (cons 1 (car Centro))
    (cons 2 (cadr Centro))
    (cons 3 (caddr Radio0))
    (cons 10
	  (distance
	    Radio0
	    (list (car Centro) (cadr Centro) (caddr Radio0))
	  )
    )
    (cons 20 Radio1)
    (cons 120 (angle Centro Radio0))
  )
)



;;; FUNCION PRINCIPAL ESPIRAL
(defun C:ESPiral (/ ESP_Puntos SuavOrig Capa EstadoCapa EstadoEcho)
  (defun *error* (msg)
    (defun *error* (msg)
      (princ "error: ")
      (princ msg)
      (princ)
    )
  )
  (alert
    (strcat "Comando: ESPIRAL" "\nAutor: ANDRES TESTI - 2000" "\nAplicaciones AutoLISP" "\nContacto: andytesti@hotmail.com")
  )
  (if (setq ESP_Puntos (ESP:TomaPuntos))
    (progn (ESP:Calculo ESP_Puntos)
	   (setq SuavOrig   (getvar "SPLINETYPE")
		 Capa	    (getvar "CLAYER")
		 EstadoCapa (entget (tblobjname "LAYER" Capa))
		 EstadoEcho (getvar "CMDECHO")
	   )
	   (entmod (subst (cons 70 0) (assoc 70 EstadoCapa) EstadoCapa)
	   )				; Descongela la Capa Actual
	   (setvar "SPLINETYPE" 5)	; Suavizado Cuadrático
	   (setvar "CMDECHO" 0)		; Uso silencioso del comando
	   (command (getcname "_PEDIT") (entlast) "s" "x")
					;Se suaviza la polylinea
	   (setvar "SPLINETYPE" SuavOrig)
	   (setvar "CMDECHO" EstadoEcho)
	   (entmod EstadoCapa)
    )
    (princ "*Comando Cancelado*")
  )
  (princ)
)
