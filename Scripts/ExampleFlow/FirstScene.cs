using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using System;

public class FirstScene : MonoBehaviour, ITrackableStateHandler
{
    private List<IAnimationExecutor> animations;
    private float timeElapsed = 0.0f;
    private bool showAnimations = false;

    void Start()
    {
        animations = new List<IAnimationExecutor>();

        var mapGameObject = transform.Find("Map").gameObject;
        var wolfCupcakeObject = transform.Find("WolfCupcake").gameObject;

        var mapDisappear = new AnimatorParameterSetter<bool>(mapGameObject.GetComponent<Animator>(), "IsDisappear", true, 9.0f, null);
        var cupcakeAppear = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsShowCupcake", true, 10.0f, null);
        var showCupcakeRecipe = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsShowRecipe", true, 11.0f, null);
        var wolfAppear = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsShowWolf", true, 17.0f, null);
        var showWolfDna = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsShowNucleotides", true, 18.0f, null);
        var color = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsColor", true, 28.5f, null);
        var showDna = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsShowDna", true, 39.8f, null);
        var geneAppear = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsGeneAppear", true, 51.8f, null);
        var noColor = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsNoColor", true, 46.5f, null);
        var showIron = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsIron", true, 73.0f, null);
        var showBananaCupcake = new AnimatorParameterSetter<bool>(wolfCupcakeObject.GetComponent<Animator>(), "IsBananaCake", true, 82.5f, null);

        animations.Add(mapDisappear);

        animations.Add(cupcakeAppear);
        animations.Add(showCupcakeRecipe);
        animations.Add(wolfAppear);
        animations.Add(showWolfDna);
        animations.Add(color);
        animations.Add(showDna);
        //animations.Add(geneAppear);
        animations.Add(noColor);
        animations.Add(showIron);
        animations.Add(showBananaCupcake);
    }

    void Update()
    {
        if (!showAnimations)
            return;
        timeElapsed += Time.deltaTime;
        /*
        if (secondsCounter > Math.Truncate(timeElapsed))
            return;
        secondsCounter++;
        */
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
