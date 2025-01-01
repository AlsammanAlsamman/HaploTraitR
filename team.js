document.addEventListener("DOMContentLoaded", function() {
  const teamContainer = document.getElementById("team-container");

  const teamMembers = [
    {
      name: "Dr. Zakaria Kehel",
      affiliation: "genetic resources scientist and senior biometrician",
      img: "man/figures/zk.jpg"
    },
    {
      name: "Alsamman Alsamman",
      affiliation: "ICARDA Bioinformatics consultant",
      img: "man/figures/ama.png"
    },
    {
      name: "Sara El-Gebali",
      affiliation: "Your Affiliation Here",
      img: "man/figures/sarah.jpg"
    }
  ];

  teamMembers.forEach(member => {
    const memberDiv = document.createElement("div");
    memberDiv.className = "team-member";
    memberDiv.style.flex = "1 1 30%";
    memberDiv.style.textAlign = "center";
    memberDiv.style.margin = "10px";

    const img = document.createElement("img");
    img.src = member.img;
    img.alt = member.name;
    img.width = 150;
    img.height = 150;
    img.style.borderRadius = "50%";

    const name = document.createElement("h3");
    name.textContent = member.name;

    const affiliation = document.createElement("p");
    affiliation.textContent = member.affiliation;

    memberDiv.appendChild(img);
    memberDiv.appendChild(name);
    memberDiv.appendChild(affiliation);

    teamContainer.appendChild(memberDiv);
  });
});
