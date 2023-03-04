using System;
using UnityEngine;
using System.Collections.Generic;
public class ColorHeatMap
{
    public ColorHeatMap()
    {
        initColorsBlocks();
    }
    private void initColorsBlocks()
    {
        ColorsOfMap.AddRange(new Color[]{
            new Color(0, 0, 0,1f),//Black
            new Color(0, 0, 1f,1f),//Blue
            new Color(0, 1f, 1f,1f),//Cyan
            new Color(0, 1f, 0,1f),//Green
            new Color(1f, 1f, 0,1f),//Yellow
            new Color(1f, 0, 0,1f),//Red
            new Color(1f, 1f, 1f,1f)// White
        });
    }
    public Color GetColorForValue(float val, float maxVal)
    {
        float valPerc = val / maxVal;// value%
        float colorPerc = 1f / (ColorsOfMap.Count-1);// % of each block of color. the last is the "100% Color"
        float blockOfColor = valPerc / colorPerc;// the integer part repersents how many block to skip
        int blockIdx = (int)Math.Truncate(blockOfColor);// Idx of 
        float valPercResidual = valPerc - (blockIdx*colorPerc);//remove the part represented of block 
        float percOfColor = valPercResidual / colorPerc;// % of color of this block that will be filled
        Color cTarget = ColorsOfMap[blockIdx];
        float deltaR = 0;
        float deltaG = 0;
        float deltaB = 0;
        if(blockIdx + 1 != 7)
        {
            Color cNext = cNext = ColorsOfMap[blockIdx + 1]; 
            deltaR =cNext.r - cTarget.r;
            deltaG =cNext.g - cTarget.g;
            deltaB =cNext.b - cTarget.b;
        }
        
        var R = cTarget.r + (deltaR * percOfColor);
        var G = cTarget.g + (deltaG * percOfColor);
        var B = cTarget.b + (deltaB * percOfColor);

        return new Color(R, G, B,1f);
    }
    public List<Color> ColorsOfMap = new List<Color>();
}