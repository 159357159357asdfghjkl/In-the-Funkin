var camNotes:FlxCamera;
function onCreatePost(){
    camNotes = new FlxCamera();
    camNotes.height = 1300;
    camNotes.bgColor = 0x00000000;
    
    FlxG.cameras.remove(game.camHUD, false);
    FlxG.cameras.remove(game.camOther, false);

    FlxG.cameras.add(camNotes, false);
    FlxG.cameras.add(game.camHUD, false);
    FlxG.cameras.add(game.camOther, false);
    game.notes.cameras = [camNotes];
    game.strumLineNotes.cameras = [camNotes];
    camNotes.angle = game.camHUD.angle;
    camNotes.zoom = game.camHUD.zoom;
    camNotes.x = game.camHUD.x;
    camNotes.y = game.camHUD.y;
}
function onUpdate(elapsed){
    if(game.camZooming)
        camNotes.zoom = game.camHUD.zoom;
}
//some basic helpers
createGlobalCallback('setCamNotes',function(x:Float,y:Float,a:Float){
camNotes.x=x; camNotes.y=y; camNotes.angle=a;
});
createGlobalCallback('round',function(v:Float):Int{
return Math.round(v);
});
createGlobalCallback('abs',function(v:Float):Float{
return Math.abs(v);
});