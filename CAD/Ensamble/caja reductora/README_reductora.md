# Reductora planetaria compacta NEMA 17 - 5:1

Diseño preliminar fabricable para usar la misma caja en el eje de base y en el eje de brazo del manipulador. La salida es coaxial con el motor, de modo que la caja puede incorporarse sin desplazar el eje geométrico de cada articulación.

## Especificación

| Parámetro | Valor |
|---|---:|
| Relación exacta | 5:1 |
| Configuración | Sol de entrada, corona fija, portaplanetas de salida |
| Engranajes | Sol 12T, 3 planetas 18T, corona 48T |
| Módulo / presión | 1,0 / 20 grados |
| Ancho dentado | 10 mm |
| Diámetro exterior | 64 mm |
| Longitud caja sin motor ni eje saliente | 26 mm |
| Entrada | Eje NEMA 17 tipo D, nominal 5 mm |
| Interfaz del motor | Resalte 22 mm y patrón 31 x 31 mm, 4 tornillos M3 |
| Salida | Eje 8 mm y rodamiento 608 (8 x 22 x 7 mm) |
| Cierre | 4 tornillos M4 en círculo de 57 mm |

La relación se obtiene de:

`i = 1 + Z_corona/Z_sol = 1 + 48/12 = 5`

## Criterio de montaje en el robot

- En la base, fijar la carcasa al NEMA 17 y unir el eje de salida de 8 mm al plato giratorio. El plato debe poseer un segundo rodamiento, separado del 608 de la tapa, para que el peso y el momento de vuelco no sean soportados por los engranajes.
- En el brazo, usar la misma caja, pero fijar el eje de salida al pivote del brazo mediante abrazadera partida o cubo con chaveta/plano. Añadir otro rodamiento del lado opuesto del brazo.
- El tercer NEMA 17 puede montarse en una placa con ranuras paralelas para tensar la correa. La polea debe sujetarse al plano del eje de 5 mm mediante dos prisioneros M3, uno de ellos enfrentado al plano.

## Componentes no impresos por caja

| Cantidad | Componente |
|---:|---|
| 1 | Rodamiento 608-2RS o 608-ZZ |
| 3 | Tornillos de hombro o pernos M5 para ejes planetarios |
| 6 | Arandelas finas M5 de baja fricción |
| 4 | Tornillos M3 para fijar al NEMA 17; longitud según motor |
| 4 | Tornillos M4 para cerrar la tapa |
| 1 | Prisionero M3 para el engranaje solar |
| 1 | Eje de acero calibrado de 8 mm para la salida definitiva |

## Fabricación recomendada

- Prototipo: PETG, 5 perímetros, 40-60 % de relleno, capas de 0,15-0,20 mm. Imprimir engranajes acostados.
- Versión de trabajo: PA-CF/POM para los engranajes o engranajes mecanizados. No usar un eje de salida impreso para cargar el brazo; el sólido integrado del STEP sirve para comprobar el conjunto y debe sustituirse por acero de 8 mm.
- Lubricar moderadamente con grasa compatible con el polímero.
- La holgura tangencial inicial es 0,18 mm. Hacer primero una probeta de sol-planeta; según la calibración de la impresora puede convenir 0,12-0,25 mm.

## Par esperado y límites

Con un NEMA 17 de 0,45 N m y una eficiencia conservadora del 70-80 %, el par de salida esperado es aproximadamente 1,6-1,8 N m. El valor real depende mucho de la velocidad, el driver, la corriente y la calidad de los dientes. Esta caja no debe absorber cargas radiales o momentos del brazo: esas cargas deben cerrarse a través de dos rodamientos del propio eje de la articulación.

Antes de fabricar la versión final hay que confirmar el par del motor, masa y longitud de cada eslabón, distancia disponible alrededor de cada articulación y diámetro/interfaz deseados en los ejes de salida.

## Regeneración

Ejecutar `python3 cad/reductora_nema17_5a1.py`. Los parámetros principales están agrupados al comienzo del archivo para modificar módulo, holgura y ancho.
