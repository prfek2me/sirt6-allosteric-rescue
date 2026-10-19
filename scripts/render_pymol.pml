reinitialize
set connect_mode, 1
load results/analysis/SIRT6_Holo_MDL800_basin_rep.pdb, sirt6_holo

bg_color white
hide all

show cartoon, polymer
color gray80, polymer

show sticks, resn MDL
util.cbac resn MDL
set stick_radius, 0.2, resn MDL

select pocket, byres polymer within 5.0 of resn MDL
show sticks, pocket
util.cbag pocket
disable pocket

distance polar_conts, resn MDL, pocket, 3.2, mode=2
hide labels, polar_conts
set dash_color, red, polar_conts
set dash_radius, 0.05, polar_conts

show sticks, resn CM1+CM2+CM3+CM4
util.cbay resn CM1+CM2+CM3+CM4

show spheres, resn ZN1
color purple, resn ZN1
set sphere_scale, 0.5, resn ZN1

distance zn_coord, resn ZN1, (resn CM1+CM2+CM3+CM4 and name SG), 2.6
hide labels, zn_coord
set dash_color, purple, zn_coord
set dash_radius, 0.06, zn_coord

zoom resn MDL or resn ZN1
