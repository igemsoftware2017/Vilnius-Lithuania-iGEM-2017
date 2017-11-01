using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.EventSystems;

public interface ITrackableStateHandler : IEventSystemHandler
{
    void OnTrackableFound(GameObject gameObject);
    void OnTrackableLost(GameObject gameObject);
}
