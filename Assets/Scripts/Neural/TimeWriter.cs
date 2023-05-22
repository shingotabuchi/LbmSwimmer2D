using System.IO;
using System.Collections;
using System.Collections.Generic;

public class TimeWriter
{
    public static string times = "";
    public static string thetas = "";
    public static void WriteTimes()
    {
        File.WriteAllText(@"pynotes/texts/grav" + GravParticle.GravScale.ToString("0.0") +"_drag" + GravParticle.DragScale.ToString() + "_initvel" +GravParticle.InitVelScale.ToString() + ".txt", times);
        times = "";

        File.WriteAllText(@"pynotes/initThetas/grav" + GravParticle.GravScale.ToString("0.0") +"_drag" + GravParticle.DragScale.ToString() + "_initvel" +GravParticle.InitVelScale.ToString() + ".txt", thetas);
        thetas = "";
    }
}