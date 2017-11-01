using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

public class ThirdScene : MonoBehaviour, ITrackableStateHandler
{
    private List<IAnimationExecutor> animations;
    private float timeElapsed = 0.0f;
    private bool showAnimations = false;

    void Start()
    {
        animations = new List<IAnimationExecutor>();

        var division = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsDivide", true, 8.0f, null);
        var movement = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsMovement", true, 11.0f, null);
        var platform = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsPlatform", true, 20.0f, null);
        var combination = new AnimatorParameterSetter<bool>(transform.GetComponent<Animator>(), "IsCombination", true, 25.0f, null);

        animations.Add(division);
        animations.Add(movement);
        animations.Add(platform);
        animations.Add(combination);
    }

    void Update()
    {
        if (!showAnimations)
            return;
        timeElapsed += Time.deltaTime;

        var animationsToUpdate = animations.Where(animation => 
            animation.GetStartTime() < timeElapsed &&
            (animation.GetEndTime() == null || animation.GetEndTime() > timeElapsed));

        if (animationsToUpdate.Count() == 0)
            return;

        animationsToUpdate.ToList().ForEach(animation => animation.Execute());

        var animationsToRemove = animations
            .Where(animation => animation.GetStartTime() < timeElapsed &&
            (animation.GetEndTime() == null || animation.GetEndTime() < timeElapsed));

        animationsToRemove.ToList().ForEach(animation => animations.Remove(animation));
    }


    public void OnTrackableFound(GameObject gameObject)
    {
        showAnimations = true;
    }

    public void OnTrackableLost(GameObject gameObject)
    {
        showAnimations = false;
    }

}
