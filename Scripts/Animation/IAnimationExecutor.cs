using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public interface IAnimationExecutor {

    void Execute();
    float GetStartTime();
    float? GetEndTime();

}
