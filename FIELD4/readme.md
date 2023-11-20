# **FIELD4**


FIELD 4 is an affordable, 4 channel open-source ambisonic microphone that could serve as an easyintroduction to mmersive audio. 




**PRACTICAL CONSTRUCTION OF THE MICROPHONE**

Building a microphone consists of three parts. The first is the actual physical construction of the microphone. The second part is the calibration of the microphone. The third part is the actual use of the microphone in the field and consequent processing of the audio captured by the microphone in a Digital Audio Workstation. While the first two parts are described here, the third part is beyond the scope of this github readme.There are plenty of excellent resources on the internet that will help one do this post-processing and conversion of A format into B format.






![Fig 1](https://github.com/theisro/UMA/assets/28617707/2fa7f6cc-91c5-49ea-9889-3382e805af23)






Figure 1: Parts needed to construct the ambisonic microphone.1. Aluminium Tube 2. Aluminium Plate 3.Template for drilling holes 4. 12 pin DIN connector 5. 14mm Electret capsules(x4) 6. Aluminium Mesh 7. Glue 8. Wires.


    1. Aluminium Tube - 100mm Length , 38mm Diameter , 34 mm inner diameter. This is the main housing container that houses the 3D printed tetrahedral array and the microphones mounted to them along with the electronics.

    2. Aluminium Plate - 10mm Thick, 40 X 40mm (Square) (Also cut at the vendor from a larger plate). This is the base component of the whole device that covers the aluminium tube housing in the bottom , acts as a protector from the elements and this is where the hole for the din is drilled. 

    3. 3D Printed Template for 24 holes. This template printed from a 3D printer slots into the main aluminium housing. It acts as a stencil to align the holes onto the housing 

    4. 12 pin DIN Connector. This is the connector that is at the base of the microphone (Mounted on the base plate). It transfers audio information that's getting inputted into the microphone. It interfaces the microphone either with recorders or audio interfaces. 

    5. 14mm Microphone capsules - - The microphone capsules used to build low cost ambisonics come in various diameters. We chose to import 14mm diameter electret capsules and use it in our build. These capsule’s cost $4 - $6 on average. You need to use 4 of these capsules to build a 4ch ambisonic microphone.

    6. Aluminium Mesh - 10.5cm X 5.5cm (Cut from a larger mesh sheet). This mesh is housed inside the aluminium housing wedged between the housing and the array. It protects the microphone capsules inside from the elements such as wind, dust etc. 
    7. Araldite Glue - To glue the aluminium housing to the base plate.
    8. Wire - You have to use the thinnest gauge wire to connect the microphone capsules to the DIN connector. These wires are routed in through the 3D printed tetrahedral array.
    9. Soldering iron
	
Tools

    1. 3D Printer - Tetrahedron Array, 24 Hole Template. 
	A 3D printer is used to print the tetrahedral array that houses the capsules inside and the 24 	hole template to use as a stencil to drill holes onto the housing.
    2. Metal Drilling Machine.
	A drilling machine with a 10mm drill bit is used to drill 24 holes on the housing.
    3. Metal Lathe.
	A lathe tool is used to trim out, shape and smoothen the housing and the base plate. It has a 	holder that spins and a sharp precise blade on the other side to trim the piece that is held by 	the holder
    4. Metal Buffing machine.
	A buffing tool is used in varying grit sizes to achieve different results , a harder/rough grit 	size to smoothen out the surface of the aluminum and a softer/smoother grit size to polish 	and give a sheen to the housing and the base.


Metal Shop

Aluminum Base and Top.


Drill a 17.5mm hole into the base of the aluminium block. Insert a bolt and a nut to attach the base to the metal lathe machine. Using the machine, trim down the base from a square to a circle of 38mm. Carve the inside of the base circle for the stepped groove to attach with the main aluminium housing. Buff using a buffing machine to smoothen the edges. Screw the 8 pin DIN through the 17.5mm hole and then glue it. Figure 1 shows the top of the microphone which can be drilled in a similar manner.




![Fig 2](https://github.com/theisro/UMA/assets/28617707/e2714673-8d39-4aef-b444-a9b75703bc71)

Figure 2: Aluminium base before and after processing






![Fig 3](https://github.com/theisro/UMA/assets/28617707/ebc83816-3416-4478-937f-40f9106bdab2)




Figure 3: Aluminium top before and after machining and drilling.



Aluminium Housing

The Aluminium housing will house the main tetrahedral capsules.
Smooth out and trim the edges of the metal housing on both the sides using the metal lathe machine. Use the 3D printed template to mark the 24 holes with a sketch pen. Use a pointed chisel to mark out the centre point of each of the 24 holes, This makes the drill bit of the drilling machine to drill the holes precisely.Use a 10mm drill bit to get the exact diameter of the hole. Once the holes are punched use a file inside and outside to remove the shards of aluminium and to smoothen the surface. Use the buffing machine with a rough grit to clean and smoothen the external surface of the aluminium using. Use the buffing machine again with a smoother grit to polish and buff the aluminium housing.




![Fig 4](https://github.com/theisro/UMA/assets/28617707/3f345d20-1835-4ebb-898c-a4f39215045b)


Figure 4: Aluminium tube with the hole template before and after drilling.





Cut a piece of mesh of the size 5.5cm X 10cm from the larger roll. Roll the piece of the cut mesh so that it slots into the metal housing. Slot the rolled piece of mesh into the main aluminum housing.S



![Fig 5](https://github.com/theisro/UMA/assets/28617707/0a0749d3-2750-4d13-94f4-54baedde374f)

Figure 5: Aluminium housing with mesh


 3D Printed templates

An STL file for the 24 holes of 10 mm is required to precisely drill a holes. A .stl file for the 24 hole template.  You can reuse this to make multiple microphones.

Use a 3D printer to  print out the tetrahedral array that holds the microphone capsules. Ideally using a high quality 3d printer.We have found out NYLON 12 is a sturdy filament for printing these files

![Fig 6](https://github.com/theisro/UMA/assets/28617707/9d22c4e7-8b9f-44d4-a09c-39770cd706a6)


Figure 6: 3D printed template for holes.



![Fig 7](https://github.com/theisro/UMA/assets/28617707/32bee21b-bd2c-4bfd-a37e-7dcf64b0ddd9)




Figure 7: 3D printed template of the tetrahedral array.


Slot the 4 capsules into the 3D printed tetradhedron and then use the circuit diagram in figure 8 to wire it to the 12 pin DIN connecter which is now in the base plate.
