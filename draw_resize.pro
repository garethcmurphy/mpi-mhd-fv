PRO Draw_Resize_TLB_Resize, event

   ; Get the info structure.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Resize the draw widget.

Widget_Control, info.drawID, Draw_XSize=event.x, Draw_YSize=event.y

   ; Make the draw window the current graphics window and load the
   ; widget colors.

WSet, info.wid
TVLCT, info.r, info.g, info.b

   ; Redisplay the graphical display.

Erase
TVImage, *info.scaled, Position=[0.1, 0.1, 0.9, 0.75]
Colorbar, Range=[Min(*info.image), Max(*info.image)], Divisions=8, $
   Title=info.title, Font=1, Position=[0.1, 0.85, 0.9, 0.90], $
	Format='(G10.2)'


   ; Store the info structure.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END ;----------------------------------------------------------------



PRO Draw_Resize_Cleanup, tlb

   ; Clean up pointers.

Widget_Control, tlb, Get_UValue=info, /No_Copy
IF N_Elements(info) EQ 0 THEN RETURN
Ptr_Free, info.image
Ptr_Free, info.scaled
END ;----------------------------------------------------------------



PRO Draw_Resize, image, title

IF N_Elements(image) EQ 0 THEN image = Dist(300, 400)

   ; Scale the image data.

scaled = BytScl(image, Top=!D.Table_Size-1)

   ; Create the widgets. Make sure the TLB returns resize events.

tlb = Widget_Base(Column=1, Title='Draw Widget Resize Example', $
   TLB_Size_Events=1)
drawID = Widget_Draw(tlb, XSize=400, YSize=400)
Widget_Control, tlb, /Realize

   ; Get the window index number. Make it current.

Widget_Control, drawID, Get_Value=wid
WSet, wid

   ; Display the data.

LoadCT, 33
TVImage, scaled, Position=[0.1, 0.1, 0.9, 0.75]
Colorbar, Range=[Min(image), Max(image)], Divisions=8, $
   Charsize=3.0, $
   Title=title, Font=0, Position=[0.1, 0.85, 0.9, 0.90]

   ; Get the colors associated with the widget so you
   ; can protect them.

TVLCT, r, g, b, /Get

   ; Create an info structure and store it.

info = {image:Ptr_New(image), scaled:Ptr_New(scaled, /No_Copy), $
   wid:wid, drawID:drawID, r:r, g:g, b:b, title:title}
Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Set up XManager.

XManager, 'draw_resize', tlb, Event_Handler='Draw_Resize_TLB_Resize', $
   Cleanup='Draw_Resize_Cleanup', /No_Block
END ;----------------------------------------------------------------
