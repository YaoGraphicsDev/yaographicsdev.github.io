---
layout: default
---

<div class="home">
  <section class="home-intro">
    <p>
      I am a graphics developer focused on rendering-engine architecture,
      physically based rendering, and simulation.
    </p>

  <h2 class="post-list-heading">Selected Projects</h2> 

    <ul>
      <li>
        <a href="https://github.com/YaoGraphicsDev/rendering_rt">
          Vulkan Rendering Engine
        </a>
      </li>
      <li>
        <a href="https://github.com/YaoGraphicsDev/otcv">
          Vulkan Framework
        </a>
      </li>
      <li>
        <a href="https://github.com/YaoGraphicsDev/featherstone">
          Articulated-Body Physics Engine
        </a>
      </li>
    </ul>
  </section>

  <h2 class="post-list-heading">Posts</h2>

  <ul class="post-list">
    {% for post in site.posts reversed %}
      <li>
        <span class="post-meta">
          {{ post.date | date: "%b %-d, %Y" }}
        </span>

        <h3>
          <a class="post-link" href="{{ post.url | relative_url }}">
            {{ post.title }}
          </a>
        </h3>
      </li>
    {% endfor %}
  </ul>
</div>